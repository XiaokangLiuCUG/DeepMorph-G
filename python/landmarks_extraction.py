# Created by Xiaokang Liu on 25/05/2026 => xkliu@cug.edu.cn

import cv2
import numpy as np
import os
import csv
from glob import glob

# ========== Configuration ==========
input_folder = r"D:\new\1"
output_csv = r"all_gastropod_landmarks.csv"
valid_extensions = ('.png', '.jpg', '.jpeg', '.bmp')

# ========== Helper functions (same as before) ==========
def imread_chinese(path):
    with open(path, 'rb') as f:
        data = f.read()
    img = cv2.imdecode(np.frombuffer(data, np.uint8), cv2.IMREAD_GRAYSCALE)
    return img

def is_clockwise(points):
    area = 0
    for i in range(len(points)):
        x1, y1 = points[i]
        x2, y2 = points[(i+1) % len(points)]
        area += (x2 - x1) * (y2 + y1)
    return area < 0

def normalize_contour(contour, target_clockwise=True):
    pts = contour.squeeze().astype(np.float32)
    if len(pts) < 3:
        return pts
    if is_clockwise(pts) != target_clockwise:
        pts = pts[::-1]
    min_y = np.min(pts[:, 1])
    candidates = np.where(pts[:, 1] == min_y)[0]
    start_idx = candidates[np.argmin(pts[candidates, 0])]
    pts = np.roll(pts, -start_idx, axis=0)
    return pts

def find_highest_point(points):
    y_min = np.min(points[:, 1])
    candidates = np.where(points[:, 1] == y_min)[0]
    if len(candidates) == 1:
        return candidates[0]
    xs = points[candidates, 0]
    sorted_idx = np.argsort(xs)
    mid_idx = sorted_idx[len(sorted_idx)//2]
    return candidates[mid_idx]

def find_lowest_point(points):
    y_max = np.max(points[:, 1])
    candidates = np.where(points[:, 1] == y_max)[0]
    if len(candidates) == 1:
        return candidates[0]
    xs = points[candidates, 0]
    sorted_idx = np.argsort(xs)
    mid_idx = sorted_idx[len(sorted_idx)//2]
    return candidates[mid_idx]

def find_left_right_at_mid(points, y_mid):
    diffs = np.abs(points[:, 1] - y_mid)
    min_diff = np.min(diffs)
    candidates = np.where(diffs == min_diff)[0]
    left_idx = candidates[np.argmin(points[candidates, 0])]
    right_idx = candidates[np.argmax(points[candidates, 0])]
    return left_idx, right_idx

def get_segment(points, start_idx, end_idx):
    n = len(points)
    if start_idx <= end_idx:
        return points[start_idx:end_idx+1]
    else:
        return np.vstack([points[start_idx:], points[:end_idx+1]])

def sample_by_index(segment_pts, n_samples, exclude_endpoints=True):
    M = len(segment_pts)
    if M < 2 or n_samples <= 0:
        return np.empty((0,2))
    if exclude_endpoints:
        indices = np.linspace(0, M-1, n_samples+2)[1:-1]
    else:
        indices = np.linspace(0, M-1, n_samples)
    indices = np.round(indices).astype(int)
    # Remove duplicates (should not occur in theory, but kept for safety)
    indices = np.unique(indices)
    # If n_samples points are required but unique leaves fewer, append the last point (rare)
    if len(indices) < n_samples and not exclude_endpoints:
        indices = np.append(indices, [M-1]*(n_samples - len(indices)))
    elif len(indices) < n_samples and exclude_endpoints:
        # If interior points are insufficient, resample proportionally
        step = M / (n_samples + 1)
        indices = np.array([int(round(step * (i+1))) for i in range(n_samples)])
    return segment_pts[indices]

def process_one_image(image_path):
    """Process a single image; return (species_name, ordered_points) or (None, None) on failure"""
    # Extract species name (first underscore-separated field of the filename)
    basename = os.path.basename(image_path)
    species = basename.split('_')[0]

    # Read image
    img = imread_chinese(image_path)
    if img is None:
        print(f"  Warning: cannot read image {basename}, skipping")
        return None, None

    # Binarize and extract the largest contour
    _, binary = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY_INV)
    contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        print(f"  Warning: no contour found in {basename}, skipping")
        return None, None

    largest_contour = max(contours, key=cv2.contourArea)
    shell_pts = normalize_contour(largest_contour, target_clockwise=True)
    shell_pts = np.round(shell_pts).astype(np.int32)

    if len(shell_pts) < 4:
        print(f"  Warning: too few contour points in {basename}, skipping")
        return None, None

    # Find the four landmarks
    idx_top = find_highest_point(shell_pts)
    idx_bottom = find_lowest_point(shell_pts)
    y_min = shell_pts[idx_top, 1]
    y_max = shell_pts[idx_bottom, 1]
    y_mid = (y_min + y_max) / 2.0
    left_idx, right_idx = find_left_right_at_mid(shell_pts, y_mid)

    # Determine clockwise order
    total = len(shell_pts)
    def clockwise_order(start, a, b, n):
        return ((a - start) % n) < ((b - start) % n)
    if not clockwise_order(idx_top, left_idx, right_idx, total):
        left_idx, right_idx = right_idx, left_idx

    # Split into segments
    seg1 = get_segment(shell_pts, idx_top, left_idx)      # P1 -> P3
    seg2 = get_segment(shell_pts, left_idx, idx_bottom)   # P3 -> P2
    seg3 = get_segment(shell_pts, idx_bottom, right_idx)  # P2 -> P4
    seg4 = get_segment(shell_pts, right_idx, idx_top)     # P4 -> P1

    # Sample semi-landmarks (20 per segment)
    n_semi = 20
    semi_seg1 = sample_by_index(seg1, n_semi, exclude_endpoints=True)
    semi_seg2 = sample_by_index(seg2, n_semi, exclude_endpoints=True)
    semi_seg3 = sample_by_index(seg3, n_semi, exclude_endpoints=True)
    semi_seg4 = sample_by_index(seg4, n_semi, exclude_endpoints=True)

    # Build the ordered point list (clockwise order)
    ordered_points = []
    # P1
    ordered_points.append(shell_pts[idx_top])
    # seg1 semi-landmarks
    for pt in semi_seg1:
        ordered_points.append(pt)
    # P3
    ordered_points.append(shell_pts[left_idx])
    # seg2 semi-landmarks
    for pt in semi_seg2:
        ordered_points.append(pt)
    # P2
    ordered_points.append(shell_pts[idx_bottom])
    # seg3 semi-landmarks
    for pt in semi_seg3:
        ordered_points.append(pt)
    # P4
    ordered_points.append(shell_pts[right_idx])
    # seg4 semi-landmarks
    for pt in semi_seg4:
        ordered_points.append(pt)

    # Check the count: 4 + 80 = 84
    if len(ordered_points) != 84:
        print(f"  Warning: {basename} incorrect number of sampled points ({len(ordered_points)}), expected 84, skipping")
        return None, None

    return species, ordered_points


# ========== Main program: traverse the folder ==========
def main():
    # Get all image paths
    image_paths = []
    for ext in valid_extensions:
        image_paths.extend(glob(os.path.join(input_folder, '*' + ext)))

    print(f"Found {len(image_paths)} images, starting processing...")

    results = []  # each element is (species, points_list)
    failed = []

    for i, img_path in enumerate(image_paths, 1):
        print(f"[{i}/{len(image_paths)}] Processing: {os.path.basename(img_path)}")
        species, points = process_one_image(img_path)
        if species is not None:
            results.append((species, points))
        else:
            failed.append(os.path.basename(img_path))

    print(f"\nProcessing complete: {len(results)} succeeded, {len(failed)} failed")
    if failed:
        print("Failed list:")
        for f in failed:
            print(f"  {f}")

    # Save summary CSV
    if not results:
        print("No images were successfully processed; exiting")
        return

    # Build column names: species, then point1_x, point1_y, point2_x, point2_y, ..., point84_x, point84_y
    header = ['species']
    for i in range(1, 85):
        header.append(f'p{i}_x')
        header.append(f'p{i}_y')

    with open(output_csv, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for species, points in results:
            # Flatten the (x,y) of 84 points into a one-dimensional list
            flat = [species]
            for (x, y) in points:
                flat.append(int(x))
                flat.append(int(y))
            writer.writerow(flat)

    print(f"Summary results saved to: {os.path.abspath(output_csv)}")

if __name__ == "__main__":
    main()