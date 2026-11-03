# Created by Xiaokang Liu on 25/06/2023 => xkliu@cug.edu.cn
import cv2
import numpy as np
from pathlib import Path
import csv

def imread_gray(path):
    with open(path, 'rb') as f:
        data = np.frombuffer(f.read(), dtype=np.uint8)
        img = cv2.imdecode(data, cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise FileNotFoundError(f"Cannot read: {path}")
    return img

def get_x_mean(mask_binary):
    ys, xs = np.where(mask_binary > 0)
    if len(xs) == 0:
        return None
    return np.mean(xs)

def get_aperture_top_bottom(aperture):
    pts = np.column_stack(np.where(aperture > 0))
    if len(pts) == 0:
        return None, None
    min_y = np.min(pts[:, 0])
    max_y = np.max(pts[:, 0])
    return min_y, max_y

def get_outline_top(outline):
    pts = np.column_stack(np.where(outline > 0))
    if len(pts) == 0:
        return None
    return np.min(pts[:, 0])

def classify_sinistral_dextral(outline_path, aperture_path):
    # Read and invert (original black foreground/white background -> white foreground/black background)
    outline_raw = imread_gray(outline_path)
    aperture_raw = imread_gray(aperture_path)
    if outline_raw.shape != aperture_raw.shape:
        raise ValueError("Shape mismatch")
    outline = cv2.bitwise_not(outline_raw)
    aperture = cv2.bitwise_not(aperture_raw)
    _, outline = cv2.threshold(outline, 127, 255, cv2.THRESH_BINARY)
    _, aperture = cv2.threshold(aperture, 127, 255, cv2.THRESH_BINARY)

    # Get the top and bottom of the aperture
    ap_top, ap_bottom = get_aperture_top_bottom(aperture)
    if ap_top is None:
        raise ValueError("Aperture region is empty")
    aperture_height = ap_bottom - ap_top
    shift = aperture_height // 2
    cut_y = max(0, ap_top - shift)

    # Get the top of the shell outline
    out_top = get_outline_top(outline)
    if out_top is None:
        raise ValueError("Outline region is empty")
    if cut_y > out_top:
        cut_y = out_top

    # Crop the lower part of the shell starting from cut_y
    bottom_outline = outline[cut_y:, :]
    outline_x_mean = get_x_mean(bottom_outline)
    if outline_x_mean is None:
        raise ValueError("No valid pixels in the lower part of the shell")

    # Compute the mean x of the aperture
    aperture_x_mean = get_x_mean(aperture)
    if aperture_x_mean is None:
        raise ValueError("No valid pixels in the aperture")

    if outline_x_mean < aperture_x_mean:
        return "dextral"
    else:
        return "sinistral"

if __name__ == "__main__":
    outline_dir = Path(r"D:\GBIFgastropoda_outline")
    aperture_dir = Path(r"D:\GBIFgastropoda_aperture")
    output_csv = "GBIF_species_data.csv"
    
    results = []
    for outline_path in outline_dir.glob("*.png"):
        stem = outline_path.stem
        aperture_path = aperture_dir / f"{stem}.png"
        if not aperture_path.exists():
            # Missing aperture mask, record as an error
            species = stem.split('_')[0]  # Species name is the part before the first underscore
            results.append([species, outline_path.name, "error_missing_aperture"])
            continue

        try:
            result = classify_sinistral_dextral(str(outline_path), str(aperture_path))
            species = stem.split('_')[0]
            results.append([species, outline_path.name, result])
        except Exception as e:
            species = stem.split('_')[0]
            results.append([species, outline_path.name, f"error_{str(e)}"])

    # Write to CSV file
    with open(output_csv, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(["species name", "image name", "shell type"])
        writer.writerows(results)

    print(f"Done! Processed {len(results)} images in total, results saved to {output_csv}")