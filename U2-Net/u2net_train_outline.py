# Created by Xiaokang Liu on 08/05/2026 => xkliu@cug.edu.cn
import os
import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.utils.data import DataLoader
from torchvision import transforms
import torch.optim as optim
from torch.optim.lr_scheduler import OneCycleLR
import numpy as np
from sklearn.model_selection import train_test_split
import torch.nn.functional as F
import csv

# Custom modules (make sure these files are in the same directory or on PYTHONPATH)
from data_loader import RescaleT, RandomCrop, ToTensorLab, SalObjDataset
from model import U2NET, U2NETP

from losses.tversky import TverskyLoss, FocalTverskyLoss,JaccardLoss,ComboLoss,FocalLoss,BoundaryLoss,BoundaryDoULoss
# ------------------------------- Loss functions -------------------------------

bce_loss = nn.BCELoss(size_average=True)  # size_average=True is deprecated in the original code, kept for compatibility
tverskyloss = TverskyLoss(ignore_index=255, smooth=1.0, alpha=0.3, beta=0.7)
focaltverskyloss= FocalTverskyLoss(ignore_index=255, smooth=1.0, alpha=0.3, beta=0.7, gamma=2.0)
jaccardloss = JaccardLoss(ignore_index=255, smooth=1.0)
comboloss = ComboLoss(ignore_index=255, smooth=1.0, alpha=0.5)
focalloss = FocalLoss(ignore_index=255, gamma=2.5, alpha=0.25)
boundaryloss = BoundaryLoss(boundary_width=3)
boundarydouloss = BoundaryDoULoss(n_classes=2, smooth=1e-5, alpha_trunc=0.8)
def save_metrics_to_csv(epoch, val_loss, results, loss_name, csv_dir="."):
    """
    Save validation metrics to a CSV file named after the loss function.
    :param epoch: current epoch number
    :param val_loss: average validation loss
    :param results: metrics dictionary
    :param loss_name: loss function name (e.g., 'tversky', 'focal_tversky', 'combo')
    :param csv_dir: directory to save the CSV file, defaults to the current directory
    """
    csv_path = os.path.join(csv_dir, f"validation_metrics_{loss_name}.csv")
    row = {
        "epoch": epoch,
        "val_loss": val_loss,
        "Overall Acc": results["Overall Acc"],
        "Mean Acc": results["Mean Acc"],
        "Mean IoU": results["Mean IoU"],
        "FreqW Acc": results["FreqW Acc"],
        "Class IoU_0": results["Class IoU"][0],
        "Class IoU_1": results["Class IoU"][1],
    }
    file_exists = os.path.isfile(csv_path)
    with open(csv_path, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=row.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

def muti_bce_loss_fusion(d0, d1, d2, d3, d4, d5, d6, labels_v):
    """
    Multi-scale BCE loss + boundary loss (only for d0)
    lambda_boundary: weight of the boundary loss, recommended 0.05 ~ 0.2
    """
    loss0 = jaccardloss(d0, labels_v)
    loss1 = jaccardloss(d1, labels_v)
    loss2 = jaccardloss(d2, labels_v)
    loss3 = jaccardloss(d3, labels_v)
    loss4 = jaccardloss(d4, labels_v)
    loss5 = jaccardloss(d5, labels_v)
    loss6 = jaccardloss(d6, labels_v)
    loss = loss0 + loss1 + loss2 + loss3 + loss4 + loss5 + loss6
    
    #Boundary loss (only for the highest-resolution d0)
    loss_boundarydou = boundarydouloss(d0, labels_v)
    loss = loss + 0.5 * loss_boundarydou #lambda_boundary
    
    return loss0, loss
    
# def muti_bce_loss_fusion(d0, d1, d2, d3, d4, d5, d6, labels_v):

    # #tverskyloss,focaltverskyloss,comboloss,jaccardloss, focalloss,bce_loss
    # loss0 = jaccardloss(d0, labels_v)*0.5 + focalloss(d0, labels_v)*0.5
    # loss1 = jaccardloss(d1, labels_v)*0.5 + focalloss(d1, labels_v)*0.5
    # loss2 = jaccardloss(d2, labels_v)*0.5 + focalloss(d2, labels_v)*0.5
    # loss3 = jaccardloss(d3, labels_v)*0.5 + focalloss(d3, labels_v)*0.5
    # loss4 = jaccardloss(d4, labels_v)*0.5 + focalloss(d4, labels_v)*0.5
    # loss5 = jaccardloss(d5, labels_v)*0.5 + focalloss(d5, labels_v)*0.5
    # loss6 = jaccardloss(d6, labels_v)*0.5 + focalloss(d6, labels_v)*0.5
    # loss = loss0 + loss1 + loss2 + loss3 + loss4 + loss5 + loss6

    # return loss0, loss
# ------------------------------- Utility functions -------------------------------
def read_txt(txt_path):
    with open(txt_path, "r", encoding='utf-8') as f:
        return [line.strip() for line in f.readlines()]
# ------------------------------- Main program -------------------------------
if __name__ == '__main__':
    # ---------- Configuration parameters ----------
    model_name = 'u2net'          # can be 'u2net' or 'u2netp'
    epoch_num = 20
    batch_size_train = 12
    batch_size_val = 1            # validation batch size is usually 1 for easier metric computation
    save_freq_epoch = 10          # save the model every 10 epochs
    loss_name = "outline_jaccardloss+0.5=loss_boundary_max_lr=0.005_div_factor=5.0, final_div_factor=10"   #########################

    # Root directory (all paths in the txt files are relative to this)
    root_dir = r"E:\Python_learning\PyTorch\U-2-Net-master"

    train_image_txt = r"outline_train_image.txt"
    train_mask_txt  = r"outline_train_mask.txt"
    val_image_txt   = r"outline_val_image.txt"
    val_mask_txt    = r"outline_val_mask.txt"

    model_dir = os.path.join(os.getcwd(), 'saved_models_fossil', model_name + os.sep)
    os.makedirs(model_dir, exist_ok=True)

    # ---------- 1. Read data (no random split) ----------
    def load_paths(txt_path, root):
        with open(txt_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f if line.strip()]
        abs_paths = []
        for p in lines:
            if os.path.isabs(p):
                abs_paths.append(p)
            else:
                abs_paths.append(os.path.join(root, p))
        return abs_paths

    train_imgs   = load_paths(train_image_txt, root_dir)
    train_labels = load_paths(train_mask_txt,  root_dir)
    val_imgs     = load_paths(val_image_txt,   root_dir)
    val_labels   = load_paths(val_mask_txt,    root_dir)

    print(f"Number of training samples: {len(train_imgs)}")
    print(f"Number of validation samples: {len(val_imgs)}")

    # Check for missing files (optional)
    for img_path in train_imgs + val_imgs:
        if not os.path.exists(img_path):
            print(f"Warning: image file does not exist {img_path}")
    for msk_path in train_labels + val_labels:
        if not os.path.exists(msk_path):
            print(f"Warning: mask file does not exist {msk_path}")

    # ---------- 2. Create Dataset and DataLoader ----------
    # Training set transform
    train_transform = transforms.Compose([
        RescaleT(320),
        #RandomCrop(288),
        ToTensorLab(flag=0)
    ])
    # Validation set transform
    val_transform = transforms.Compose([
        RescaleT(320),
        ToTensorLab(flag=0)
    ])

    train_dataset = SalObjDataset(
        img_name_list=train_imgs,
        lbl_name_list=train_labels,
        transform=train_transform
    )
    val_dataset = SalObjDataset(
        img_name_list=val_imgs,
        lbl_name_list=val_labels,
        transform=val_transform
    )

    train_loader = DataLoader(train_dataset, batch_size=batch_size_train,
                              shuffle=True, num_workers=1)
    val_loader = DataLoader(val_dataset, batch_size=batch_size_val,
                            shuffle=False, num_workers=1)

    # ---------- 3. Define the model ----------
    if model_name == 'u2net':
        net = U2NET(3, 1)
    elif model_name == 'u2netp':
        net = U2NETP(3, 1)
    else:
        raise ValueError("model_name must be 'u2net' or 'u2netp'")

    if torch.cuda.is_available():
        net.cuda()
        print("Training on GPU")
    else:
        print("Training on CPU")

    # ---------- 4. Optimizer ----------
    optimizer = optim.Adam(net.parameters(), lr=0.0001, betas=(0.9, 0.999),
                           eps=1e-08, weight_decay=0)
    steps_per_epoch = len(train_loader)
    total_steps = epoch_num * steps_per_epoch   # total number of training steps
    scheduler = OneCycleLR(optimizer, max_lr=0.005, total_steps=total_steps,   # total steps  # peak learning rate (adjust based on experience)
                           pct_start=0.3, anneal_strategy='cos', # first 30% of steps used for learning rate warm-up, then decay (default 0.3) # cosine annealing
                           div_factor=5.0, final_div_factor=10)
    #                                 10                    25
    # #Create scheduler: multiply learning rate by 0.5 every 5 epochs
    # scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.8)

    # ---------- 5. Load pre-trained weights (optional) ----------
    pretrain_path = r'E:\Python_learning\PyTorch\U-2-Net-master\saved_models\u2net\u2net.pth'
    if os.path.isfile(pretrain_path):
        net.load_state_dict(torch.load(pretrain_path, map_location='cpu'))
        print("Successfully loaded pre-trained weights")
    else:
        print("Pre-trained weights not found, training from scratch")
    print(loss_name)
    print("Using dynamic learning rate")
    # ---------- 6. Training and validation loop ----------
    print("--- Start training ---")
    ite_num = 0
    running_loss = 0.0
    running_tar_loss = 0.0
    ite_num4val = 0
    best_iou_class1 = 0.0
    
    for epoch in range(epoch_num):
        net.train()
        train_loss_epoch = 0.0
        train_tar_epoch = 0.0
        num_train_batches = 0

        for i, data in enumerate(train_loader):
            ite_num += 1
            ite_num4val += 1
            num_train_batches += 1

            inputs, labels = data['image'], data['label']
            inputs = inputs.type(torch.FloatTensor)
            labels = labels.type(torch.FloatTensor)

            if torch.cuda.is_available():
                inputs_v = Variable(inputs.cuda(), requires_grad=False)
                labels_v = Variable(labels.cuda(), requires_grad=False)
            else:
                inputs_v = Variable(inputs, requires_grad=False)
                labels_v = Variable(labels, requires_grad=False)

            optimizer.zero_grad()
            d0, d1, d2, d3, d4, d5, d6 = net(inputs_v)
            loss2, loss = muti_bce_loss_fusion(d0, d1, d2, d3, d4, d5, d6, labels_v)

            loss.backward()
            optimizer.step()
            scheduler.step()
            
            running_loss += loss.item()
            running_tar_loss += loss2.item()
            train_loss_epoch += loss.item()
            train_tar_epoch += loss2.item()

            # Print training info every 40 iterations
            if ite_num % 40 == 0:
                avg_loss = running_loss / ite_num4val
                avg_tar = running_tar_loss / ite_num4val
                print(f"[epoch: {epoch+1:3d}/{epoch_num:3d}, "
                      f"batch: {(i+1)*batch_size_train:5d}/{len(train_imgs):5d}, "
                      f"ite: {ite_num:6d}] train loss: {avg_loss:.6f}, tar: {avg_tar:.6f}")

            # Free intermediate variables to save GPU memory
            del d0, d1, d2, d3, d4, d5, d6, loss2, loss

        # End of one training epoch; print average training loss
        avg_train_loss = train_loss_epoch / num_train_batches
        avg_train_tar = train_tar_epoch / num_train_batches
        print(f"\n===== Epoch {epoch+1} training finished =====")
        print(f"Average training loss: {avg_train_loss:.6f}, average tar loss: {avg_train_tar:.6f}")

        # ---------- Validation phase ----------
        # Ensure StreamSegMetrics is imported at the top of the file
        from stream_metrics import StreamSegMetrics
        
        net.eval()
        metrics = StreamSegMetrics(n_classes=2)   # create evaluator for binary classification
        val_loss_sum = 0.0
        num_val_batches = 0
        best_epoch = 0
        
        with torch.no_grad():
            for i, data in enumerate(val_loader):
                inputs, labels = data['image'], data['label']
                inputs = inputs.type(torch.FloatTensor)
                labels = labels.type(torch.FloatTensor)

                if torch.cuda.is_available():
                    inputs_v = inputs.cuda()
                    labels_v = labels.cuda()
                else:
                    inputs_v = inputs
                    labels_v = labels

                d0, d1, d2, d3, d4, d5, d6 = net(inputs_v)
                _, loss = muti_bce_loss_fusion(d0, d1, d2, d3, d4, d5, d6, labels_v)
                val_loss_sum += loss.item()

                # Convert model output to binary prediction (B,1,H,W) -> (B,H,W)
                pred = (d0 > 0.5).long().squeeze(1)
                # Handle label format
                if labels_v.dim() == 4 and labels_v.shape[1] == 1:
                    target = labels_v.squeeze(1).long()
                else:
                    target = labels_v.long()

                # Update confusion matrix sample by sample
                for lt, lp in zip(target.cpu().numpy(), pred.cpu().numpy()):
                    metrics.update([lt], [lp])

                num_val_batches += 1

        # Compute average validation loss and retrieve all evaluation metrics
        avg_val_loss = val_loss_sum / num_val_batches
        results = metrics.get_results()

        print(f"===== Validation results =====")
        print(f"Average validation loss: {avg_val_loss:.6f}")
        print(f"Overall Acc: {results['Overall Acc']:.4f}")
        print(f"Mean Acc: {results['Mean Acc']:.4f}")
        print(f"Mean IoU: {results['Mean IoU']:.4f}")
        print(f"FreqW Acc: {results['FreqW Acc']:.4f}")
        print("Class IoU:")
        for cls, iou in results['Class IoU'].items():
            print(f"  Class {cls}: {iou:.4f}")
        save_metrics_to_csv(epoch+1, avg_val_loss, results, loss_name, "./results")
        
        # ---------- Save model ----------
        if (epoch + 1) % save_freq_epoch == 0:
            save_path = os.path.join(model_dir, f"{model_name}_{loss_name}.pth")
            torch.save(net.state_dict(), save_path)
            print(f"Model saved to: {save_path}")
        
        # ---------- Save best model ----------
        class1_iou = results['Class IoU'][1]
        if class1_iou > best_iou_class1:
            best_iou_class1 = class1_iou
            best_epoch = epoch + 1
            best_model_path = os.path.join(model_dir, f"{model_name}_{loss_name}_best_class1_iou.pth")
            torch.save(net.state_dict(), best_model_path)
            print(f"✓ Best model updated: Class1 IoU = {class1_iou:.4f} (Epoch {best_epoch})")

        # Reset iteration counters (consistent with original code style)
        running_loss = 0.0
        running_tar_loss = 0.0
        ite_num4val = 0

    print("Training complete!")

'''
# Compute the number of iterations per epoch (i.e., the length of the dataloader)
steps_per_epoch = len(salobj_dataloader)  # because your dataloader variable is named salobj_dataloader
total_steps = epoch_num * steps_per_epoch   # total number of training steps

# Create the OneCycleLR scheduler
# max_lr is usually set to 5-10x the initial learning rate (or higher); e.g., if the original lr=0.0001, then max_lr=0.001 or 0.005
scheduler = OneCycleLR(
    optimizer,
    max_lr=0.001,              # peak learning rate (adjust based on experience)
    total_steps=total_steps,   # total steps
    pct_start=0.3,             # first 30% of steps used for learning rate warm-up, then decay (default 0.3)
    anneal_strategy='cos',     # cosine annealing
    div_factor=25.0,           # initial learning rate = max_lr / div_factor
    final_div_factor=1e4       # final learning rate = max_lr / (div_factor * final_div_factor)
)
Parameter descriptions (for your scenario):

max_lr: the highest learning rate you want to reach. The original fixed learning rate was 0.0001; try 0.001 or 0.005.

total_steps: total number of training steps = epoch_num x number of batches per epoch.

pct_start: the proportion of total steps used for the learning rate warm-up phase (default 0.3, i.e., the first 30% of steps increase, then the remaining 70% decrease).

div_factor: initial learning rate = max_lr / div_factor, default 25. For example, max_lr=0.001 -> initial learning rate = 0.001/25 = 4e-5.

final_div_factor: final learning rate = max_lr / (div_factor * final_div_factor), default 1e4. For example, with the above settings the final learning rate is about 4e-9 (very small, helps convergence).

'''