# Created by Xiaokang Liu on 16/11/2025 => xkliu@cug.edu.cn

#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Time : 2025/11/14
# @Author : XiaoKang_Liu

import os
import time
import random
import requests
from urllib.parse import urljoin
import csv
from lxml import etree

headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 '
                  '(KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
}

def make_session():
    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                      "(KHTML, like Gecko) Chrome/116.0 Safari/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "en-GB,en;q=0.9",
        "Referer": "https://conchology.be/",
        "Connection": "keep-alive"
    })
    return session

def fetch_with_retry(session, url, max_retry=5):
    for attempt in range(max_retry):
        resp = session.get(url, timeout=220)
        if resp.status_code == 200:
            return resp
        else:
            print(f"Attempt {attempt+1} for {url} returned {resp.status_code}")
            # If 403, likely rate-limited; back off and retry
            if resp.status_code == 403:
                sleep = random.uniform(5, 8)
                print(f"Got 403, sleeping for {sleep:.1f} seconds")
                time.sleep(sleep)
            else:
                time.sleep(random.uniform(1, 3))
    # Final attempt
    return session.get(url, timeout=220)

def extract_family_links(family_url):
    session = make_session()
    resp = fetch_with_retry(session, family_url)
    resp.raise_for_status()
    tree = etree.HTML(resp.content)
    hrefs = tree.xpath(
        '//td[@class="view" and @align="right"]/a[@class="nbsp-20 hidden-xs" and contains(@href, "show=allimages")]/@href'
    )
    return [urljoin(family_url, href) for href in hrefs]

def extract_species_links(family_links):
    session = make_session()
    species_links = set()

    for fam in family_links:
        fam = fam.replace(
            "genus=&commonNames=&locality=&species=&book=&shell_ID=&form=&view=&Genus_ID",
            "p"
        )
        for i in range(20):
            resp = fetch_with_retry(session, fam + str(i + 1))
            try:
                if resp.status_code != 200:
                    print(f"Skip {fam} because status {resp.status_code}")
                    continue
                tree = etree.HTML(resp.content)
                hrefs = tree.xpath(
                    '//div[contains(@class,"data-display")]/a[1]'
                    '[starts-with(@href, "?t=263") and contains(@href,"fullspecies=") and contains(@href,"shellID=")]/@href'
                )
                if len(hrefs) == 0:
                    break
                else:
                    print("Number of species on current page:", len(hrefs))
                    for href in hrefs:
                        species_links.add(urljoin(fam, href))
                # Delay to avoid rate limiting
                time.sleep(random.uniform(3, 5))
            except Exception as e:
                print(f"Exception while processing page {fam + str(i + 1)}: {e}")
                continue  # Continue to next page on exception
    return list(species_links)

def save_links_to_file(links, filename="family_links.txt"):
    """Save links to a file."""
    try:
        with open(filename, 'w', encoding='utf-8') as f:
            for link in links:
                f.write(link + '\n')
        print(f"Successfully saved {len(links)} links to {filename}")
    except Exception as e:
        print(f"Error saving file: {e}")

if __name__ == "__main__":
    base_family_url = "https://conchology.be/?t=261&family="
    all_family_links = []
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    for i, letter in enumerate(letters, 1):
        url = base_family_url + letter
        print(f"\n[{i}/26] Processing letter {letter}")
        family_links = extract_family_links(url)
        print(f"Extracted {len(family_links)} links from {url}")
        all_family_links.extend(family_links)
        time.sleep(random.uniform(5, 7))

    print(f"Total family links found: {len(all_family_links)}")
    print("all_family_links len:", len(all_family_links))

    species_links = extract_species_links(all_family_links)
    save_links_to_file(species_links, filename="species_links.txt")
    print("species_links len:", len(species_links))