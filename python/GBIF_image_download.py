# Created by Xiaokang Liu on 12/04/2026 => xkliu@cug.edu.cn
import os
import re
import time
import random
import string
import requests
from urllib.parse import urlparse
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
import csv

# ==================== Configuration ====================
txt_file = r"D:\selected_gastropoda_image_urls.txt"
save_root = r"D:\GBIFgastropoda_images"
log_csv = r"D:\download_log.csv"

MAX_WORKERS = 8
DELAY_BETWEEN_SUBMIT = 0.1
TIMEOUT = 300
RETRY = 4
HEADERS = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'}

# Global locks
csv_lock = Lock()
status_lock = Lock()

# ==================== Helper functions ====================
def sanitize_filename(name):
    illegal_chars = r'[<>:"/\\|?*]'
    name = re.sub(illegal_chars, '_', name).strip().strip('.')
    return name[:200] if len(name) > 200 else name

def extract_original_filename(url):
    parsed = urlparse(url)
    path = parsed.path
    if '?' in path:
        path = path.split('?')[0]
    base = os.path.basename(path)
    if '.' in base:
        name, ext = os.path.splitext(base)
        return name, ext
    else:
        return base, ''

def get_extension_from_url(url):
    path = urlparse(url).path
    ext = os.path.splitext(path)[1]
    if ext and len(ext) <= 5 and ext[0] == '.':
        return ext.lower()
    return None

def random_char():
    return random.choice(string.ascii_lowercase)

def append_log(row):
    """Thread-safe append of one row to the CSV."""
    with csv_lock:
        with open(log_csv, 'a', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            writer.writerow(row)
            f.flush()

def load_existing_log():
    """Load the existing log and return a url -> (status, img_name, species) mapping."""
    url_info = {}
    if not os.path.exists(log_csv):
        with open(log_csv, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            writer.writerow(['species_name', 'url', 'img_name', 'download_status'])
        return url_info
    with open(log_csv, 'r', newline='', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            if len(row) >= 4:
                species, url, img_name, status = row[0], row[1], row[2], row[3]
                if url not in url_info:
                    url_info[url] = (status, img_name, species)
    return url_info

def download_one(url, species):
    """Download a single URL and return (status, img_name)."""
    # Build save path
    base_name_no_ext, ext_from_url = extract_original_filename(url)
    if not base_name_no_ext:
        base_name_no_ext = "image"
    ext = ext_from_url if ext_from_url else get_extension_from_url(url)
    if not ext:
        ext = '.jpg'
    elif not ext.startswith('.'):
        ext = '.' + ext
    if base_name_no_ext.lower().endswith(ext.lower()):
        base_name_no_ext = base_name_no_ext[:-len(ext)]
    species_safe = sanitize_filename(species)
    species_dir = os.path.join(save_root, species_safe)
    os.makedirs(species_dir, exist_ok=True)
    rand_char = random_char()
    final_filename = f"{species_safe}_{base_name_no_ext}_{rand_char}{ext}"
    save_path = os.path.join(species_dir, final_filename)

    if os.path.exists(save_path):
        return 'skip_exist', final_filename

    for attempt in range(RETRY):
        try:
            resp = requests.get(url, headers=HEADERS, timeout=TIMEOUT, stream=True)
            resp.raise_for_status()
            with open(save_path, 'wb') as f:
                for chunk in resp.iter_content(8192):
                    if chunk:
                        f.write(chunk)
            return 'success', final_filename
        except Exception as e:
            if attempt < RETRY - 1:
                time.sleep(2)
            else:
                return 'failed', None
    return 'failed', None

# ==================== Main program ====================
# 1. Read the text file
print("Reading text file...")
all_tasks = []  # (species, url)
with open(txt_file, 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line or '|' not in line:
            continue
        parts = line.split('|', 1)
        if len(parts) != 2:
            continue
        species, url = parts[0].strip(), parts[1].strip()
        if url:
            all_tasks.append((species, url))
total_tasks = len(all_tasks)
print(f"Total tasks (including duplicate URLs): {total_tasks}")
random.seed(42)
# 2. Load existing log (resume support)
url_info = load_existing_log()
print(f"Recovered {len(url_info)} URL records from the log")

# 3. Build a mapping from URL to its first species for duplicate detection
first_species_for_url = {}
for species, url in all_tasks:
    if url not in first_species_for_url:
        first_species_for_url[url] = species

# 4. Process tasks
def process_task(species, url):
    first_sp = first_species_for_url[url]
    # If the URL is a duplicate and the first task already succeeded or was skipped, mark as duplicate
    if species != first_sp:
        if url in url_info:
            status, img_name, _ = url_info[url]
            if status in ('success', 'skip_exist'):
                append_log([species, url, img_name, 'duplicate'])
                return
    # Otherwise download (or handle the first occurrence of the URL)
    status, img_name = download_one(url, species)
    # Update the in-memory url_info to avoid duplicate downloads
    with status_lock:
        if url not in url_info:
            url_info[url] = (status, img_name, species)
    append_log([species, url, img_name if img_name else '', status])

# Multi-threaded execution
print(f"Starting processing with {MAX_WORKERS} threads")
completed = 0
completed_lock = Lock()

def callback(future):
    global completed
    with completed_lock:
        completed += 1
        if completed % 1000 == 0:
            print(f"Progress: {completed}/{total_tasks} tasks completed")

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = []
    for species, url in all_tasks:
        future = executor.submit(process_task, species, url)
        future.add_done_callback(callback)
        futures.append(future)
        time.sleep(DELAY_BETWEEN_SUBMIT)
    # Wait for all tasks to finish
    for future in futures:
        future.result()

print("All tasks completed!")

# 5. Final statistics
print("\nComputing final statistics...")
# Read the full log file for statistics
df_log = []
with open(log_csv, 'r', newline='', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    header = next(reader)
    for row in reader:
        if len(row) >= 4:
            df_log.append(row)

success_cnt = sum(1 for r in df_log if r[3] == 'success')
duplicate_cnt = sum(1 for r in df_log if r[3] == 'duplicate')
skip_exist_cnt = sum(1 for r in df_log if r[3] == 'skip_exist')
failed_cnt = sum(1 for r in df_log if r[3] == 'failed')

print(f"\nFinal statistics:")
print(f"  Total records: {len(df_log)}")
print(f"  Successfully downloaded: {success_cnt}")
print(f"  Duplicate URLs skipped: {duplicate_cnt}")
print(f"  Existing files skipped: {skip_exist_cnt}")
print(f"  Download failed: {failed_cnt}")