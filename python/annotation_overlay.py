# Created by Xiaokang Liu on 13/04/2026 => xkliu@cug.edu.cn
# -*- coding: utf-8 -*-
import cv2
import numpy as np
from pathlib import Path
import os

def imread_chinese(path):
    with open(path, 'rb') as f:
        data = np.frombuffer(f.read(), dtype=np.uint8)
        img = cv2.imdecode(data, cv2.IMREAD_COLOR)
    if img is None:
        raise ValueError(f"无法读取: {path}")
    return img

def imread_grayscale_chinese(path):
    with open(path, 'rb') as f:
        data = np.frombuffer(f.read(), dtype=np.uint8)
        img = cv2.imdecode(data, cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise ValueError(f"无法读取: {path}")
    return img

def overlay_mask(image_path, mask_path, output_path, alpha=0.5, color=(0, 0, 255), invert=False):
    img = imread_chinese(image_path)
    mask = imread_grayscale_chinese(mask_path)
    
    if mask.shape[:2] != img.shape[:2]:
        mask = cv2.resize(mask, (img.shape[1], img.shape[0]), interpolation=cv2.INTER_NEAREST)
    
    _, binary_mask = cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY)
    
    if invert:
        binary_mask = cv2.bitwise_not(binary_mask)
    
    mask_foreground = binary_mask > 0
    overlay = np.zeros_like(img)
    overlay[:] = color
    result = img.copy()
    result[mask_foreground] = (1 - alpha) * img[mask_foreground] + alpha * overlay[mask_foreground]
    result = np.clip(result, 0, 255).astype(np.uint8)
    
    success, buf = cv2.imencode('.png', result)
    if success:
        with open(output_path, 'wb') as f:
            f.write(buf)
        print(f"已保存: {output_path}")
    else:
        raise IOError(f"保存失败: {output_path}")

def batch_overlay(image_dir, mask_dir, output_dir, alpha=0.5, color=(0,0,255), mask_suffix="_mask.png", invert=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    image_extensions = {'.jpg', '.jpeg', '.png', '.bmp', '.tif', '.JPG', '.JPEG', '.PNG', '.BMP', '.TIF'}
    image_paths = [p for p in Path(image_dir).iterdir() if p.suffix in image_extensions]
    
    for img_path in image_paths:
        stem = img_path.stem
        possible_names = [stem + mask_suffix, stem + '.png', stem + '.jpg', stem + '_mask.png', stem + '_mask.jpg', stem + '_mask' + img_path.suffix]
        mask_path = None
        for name in possible_names:
            candidate = Path(mask_dir) / name
            if candidate.exists():
                mask_path = candidate
                break
        if mask_path is None:
            print(f"跳过 {img_path.name}，未找到mask")
            continue
        output_path = output_dir / (stem + "_overlay.png")
        overlay_mask(str(img_path), str(mask_path), str(output_path), alpha, color, invert)

if __name__ == "__main__":
    image_dir = r"D:\项目\DeepMorph-G\single"
    mask_dir = r"E:\Python_learning\PyTorch\U-2-Net-master\test_data\gastropods"
    output_dir = r"D:\项目\DeepMorph-G\overlay_results"
    
    # 如果叠加区域反了，把 invert=False 改为 invert=True
    batch_overlay(image_dir, mask_dir, output_dir, alpha=0.5, color=(0,0,255), invert=True)