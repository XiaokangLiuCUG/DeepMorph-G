import torch
import torch.nn as nn
from losses.utils import flatten
import numpy as np
import torch.nn.functional as F

class TverskyLoss(nn.Module):
    def __init__(self, ignore_index=255, smooth=1.0, alpha=0.5, beta=0.5):
        super(TverskyLoss, self).__init__()
        self.ignore_index = ignore_index
        self.smooth = smooth
        self.alpha = alpha
        self.beta = beta

    def forward(self, input, target):
        # input: (B,1,H,W) 概率图 (已经 sigmoid)
        # target: (B,H,W) 长整型，0=背景，1=壳口
        input, target = flatten(input, target, self.ignore_index)  # input shape (N, 1), target shape (N,)
        
        # 对于二分类，将单通道概率转为两个通道的概率： [p, 1-p]
        prob_foreground = input[:, 0]          # (N,)
        prob_background = 1 - prob_foreground  # (N,)
        prob = torch.stack([prob_background, prob_foreground], dim=1)  # (N, 2)
        
        num_classes = 2
        losses = []
        for c in range(num_classes):
            target_c = (target == c).float()   # (N,)
            input_c = prob[:, c]               # (N,)
            tp = (input_c * target_c).sum()
            fp = ((1 - target_c) * input_c).sum()
            fn = (target_c * (1 - input_c)).sum()
            tversky = (tp + self.smooth) / (tp + self.alpha * fp + self.beta * fn + self.smooth)
            losses.append(1 - tversky)
        
        return torch.stack(losses).mean()


class FocalTverskyLoss(nn.Module):
    def __init__(self, ignore_index=255, smooth=1.0, alpha=0.3, beta=0.7, gamma=2.0):
        """
        Focal Tversky Loss for binary segmentation
        :param ignore_index: label value to ignore
        :param smooth: smoothing factor
        :param alpha: weight for false positives (FP), recommended 0.3
        :param beta:  weight for false negatives (FN), recommended 0.7 (increase penalty for missed detections)
        :param gamma: focal exponent, recommended 2.0
        """
        super(FocalTverskyLoss, self).__init__()
        self.ignore_index = ignore_index
        self.smooth = smooth
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

    def forward(self, input, target):
        """
        :param input: model output, shape (B,1,H,W), already passed through sigmoid, values in (0,1) as foreground probability
        :param target: ground truth labels, shape (B,H,W), values 0 (background) or 1 (foreground), ignore_index for ignored pixels
        """
        # 1. Flatten and filter ignore_index
        input, target = flatten(input, target, self.ignore_index)   # input: (N, 1), target: (N,)
        
        # 2. Construct two-channel probabilities: background probability = 1 - foreground probability, foreground probability = input[:,0]
        prob_foreground = input[:, 0]               # (N,)
        prob_background = 1 - prob_foreground       # (N,)
        prob = torch.stack([prob_background, prob_foreground], dim=1)  # (N, 2)
        
        # 3. Compute Focal Tversky Loss for each class
        num_classes = 2
        losses = []
        for c in range(num_classes):
            target_c = (target == c).float()        # (N,)
            input_c = prob[:, c]                   # (N,)
            
            tp = (input_c * target_c).sum()
            fp = ((1 - target_c) * input_c).sum()
            fn = (target_c * (1 - input_c)).sum()
            
            tversky = (tp + self.smooth) / (tp + self.alpha * fp + self.beta * fn + self.smooth)
            focal_tversky = (1 - tversky) ** self.gamma
            losses.append(focal_tversky)
        
        # 4. Return average loss over all classes
        return torch.stack(losses).mean()


class DiceLoss(nn.Module):
    def __init__(self, ignore_index=255, smooth=1.0):
        """
        Dice Loss for binary segmentation
        :param ignore_index: label value to ignore
        :param smooth: smoothing factor
        """
        super(DiceLoss, self).__init__()
        self.ignore_index = ignore_index
        self.smooth = smooth

    def forward(self, input, target):
        """
        :param input: model output, shape (B,1,H,W), already passed through sigmoid, values in (0,1) as foreground probability
        :param target: ground truth labels, shape (B,H,W), values 0 (background) or 1 (foreground), ignore_index for ignored pixels
        """
        # 1. Flatten and filter ignore_index
        input, target = flatten(input, target, self.ignore_index)   # input: (N, 1), target: (N,)

        # 2. Construct two-channel probabilities: background prob = 1 - foreground prob, foreground prob = input[:,0]
        prob_foreground = input[:, 0]               # (N,)
        prob_background = 1 - prob_foreground       # (N,)
        prob = torch.stack([prob_background, prob_foreground], dim=1)  # (N, 2)

        # 3. Compute Dice Loss for each class
        num_classes = 2
        losses = []
        for c in range(num_classes):
            target_c = (target == c).float()        # (N,)
            input_c = prob[:, c]                   # (N,)

            intersection = (input_c * target_c).sum()
            dice = (2. * intersection + self.smooth) / (input_c.sum() + target_c.sum() + self.smooth)

            losses.append(1 - dice)

        # 4. Return average loss over all classes
        return torch.stack(losses).mean()


class JaccardLoss(nn.Module):
    def __init__(self, ignore_index=255, smooth=1.0):
        """
        Jaccard (IoU) Loss for binary segmentation
        :param ignore_index: label value to ignore
        :param smooth: smoothing factor
        """
        super(JaccardLoss, self).__init__()
        self.ignore_index = ignore_index
        self.smooth = smooth

    def forward(self, input, target):
        """
        :param input: model output, shape (B,1,H,W), already passed through sigmoid, values in (0,1) as foreground probability
        :param target: ground truth labels, shape (B,H,W), values 0 (background) or 1 (foreground), ignore_index for ignored pixels
        """
        # 1. Flatten and filter ignore_index
        input, target = flatten(input, target, self.ignore_index)   # input: (N, 1), target: (N,)

        # 2. Construct two-channel probabilities: background prob = 1 - foreground prob, foreground prob = input[:,0]
        prob_foreground = input[:, 0]               # (N,)
        prob_background = 1 - prob_foreground       # (N,)
        prob = torch.stack([prob_background, prob_foreground], dim=1)  # (N, 2)

        # 3. Compute Jaccard Loss for each class
        num_classes = 2
        losses = []
        for c in range(num_classes):
            target_c = (target == c).float()        # (N,)
            input_c = prob[:, c]                   # (N,)

            intersection = (input_c * target_c).sum()
            total = (input_c + target_c).sum()
            union = total - intersection
            iou = (intersection + self.smooth) / (union + self.smooth)

            losses.append(1 - iou)

        # 4. Return average loss over all classes
        return torch.stack(losses).mean()


class ComboLoss(nn.Module):
    def __init__(self, ignore_index=255, smooth=1.0, alpha=0.5):
        """
        Combo Loss = alpha * Dice Loss + (1 - alpha) * Weighted Cross-Entropy Loss
        Adapted for binary segmentation with single-channel sigmoid output.
        :param ignore_index: label value to ignore
        :param smooth: smoothing factor for Dice
        :param alpha: weight for Dice loss (1-alpha for CE loss)
        """
        super(ComboLoss, self).__init__()
        self.ignore_index = ignore_index
        self.smooth = smooth
        self.alpha = alpha

    def forward(self, input, target):
        """
        :param input: model output, shape (B,1,H,W), already passed through sigmoid, values in (0,1) as foreground probability
        :param target: ground truth labels, shape (B,H,W), values 0 (background) or 1 (foreground), ignore_index for ignored pixels
        """
        # 1. Flatten and filter ignore_index (for Dice part)
        input_flat, target_flat = flatten(input, target, self.ignore_index)  # input_flat: (N,1), target_flat: (N,)

        # 2. Dice loss part (using two-class probability)
        prob_foreground = input_flat[:, 0]          # (N,)
        prob_background = 1 - prob_foreground
        prob = torch.stack([prob_background, prob_foreground], dim=1)  # (N,2)

        num_classes = 2
        dice_losses = []
        for c in range(num_classes):
            target_c = (target_flat == c).float()
            input_c = prob[:, c]
            intersection = (input_c * target_c).sum()
            dice = (2. * intersection + self.smooth) / (input_c.sum() + target_c.sum() + self.smooth)
            dice_losses.append(1 - dice)
        dice_loss = torch.stack(dice_losses).mean()

        # 3. Weighted Cross-Entropy loss part
        # Convert single-channel sigmoid output to two-channel logits (inverse sigmoid) for CE
        # Since input is already sigmoid, we compute logits = log(input/(1-input)) but to avoid numerical issues,
        # we can directly use BCE with logits by taking inverse. Simpler: use BCEWithLogitsLoss with logits = inverse_sigmoid.
        # However, we can also treat it as a standard CE on two-class probabilities using NLL loss.
        # To keep it simple and stable, we compute binary cross-entropy on the foreground channel.
        # We'll use the flattened input_flat (foreground prob) and target_flat (0/1) with a weight for positive class.
        
        # Compute class weights (inverse frequency) to handle imbalance (optional)
        # For shell aperture, foreground is small -> give higher weight to foreground.
        # Let's compute weight for foreground = (1 - fg_ratio) / fg_ratio, but we can set manually.
        # Here we provide a simple version: weight = [1.0, w_pos] where w_pos > 1.
        # For demonstration, we set w_pos = 1.0 (no weighting) but you can adjust.
        pos_weight = torch.tensor([1.0], device=input.device)  # Increase for more recall
        
        # BCE loss (without reduction, then mask ignore_index)
        bce_loss = F.binary_cross_entropy(input_flat[:, 0], target_flat.float(), reduction='none')
        # Mask ignore_index
        mask = (target_flat != self.ignore_index)
        bce_loss = (bce_loss * mask).sum() / (mask.sum() + 1e-8)
        
        # Alternative: using weighted CE with two-class probabilities (more standard)
        # But above is simpler and equivalent for binary.
        # Let's use the BCE version which directly works with sigmoid output.
        
        # 4. Combo loss
        combo = self.alpha * dice_loss + (1 - self.alpha) * bce_loss
        return combo


class FocalLoss(nn.Module):
    def __init__(self, ignore_index=255, gamma=2.0, alpha=0.25):
        """
        Focal Loss for binary segmentation
        :param ignore_index: label value to ignore
        :param gamma: focusing parameter (default 2.0)
        :param alpha: balancing factor for positive class (foreground), default 0.25
                      if None, no class balancing.
        """
        super(FocalLoss, self).__init__()
        self.ignore_index = ignore_index
        self.gamma = gamma
        self.alpha = alpha

    def forward(self, input, target):
        # input: (B,1,H,W) probability map (already sigmoid)
        # target: (B,H,W) long integer, 0=background, 1=shell aperture
        input, target = flatten(input, target, self.ignore_index)  # input shape (N,1), target shape (N,)

        # For binary segmentation, convert single-channel probability to two channels: [background_prob, foreground_prob]
        prob_foreground = input[:, 0]          # (N,)
        prob_background = 1 - prob_foreground  # (N,)
        prob = torch.stack([prob_background, prob_foreground], dim=1)  # (N, 2)

        num_classes = 2
        losses = []
        for c in range(num_classes):
            target_c = (target == c).float()   # (N,)
            input_c = prob[:, c]               # (N,)

            # p_t: predicted probability for the true class
            p_t = torch.where(target_c == 1, input_c, 1 - input_c)

            # Standard binary cross entropy term
            ce_loss = -torch.log(p_t + 1e-8)   # avoid log(0)

            # Focal weight (1 - p_t)^gamma
            focal_weight = (1 - p_t) ** self.gamma

            # Class balancing weight (alpha for foreground, 1-alpha for background)
            if self.alpha is not None:
                alpha_factor = torch.where(target_c == 1, self.alpha, 1 - self.alpha)
                focal_loss = alpha_factor * focal_weight * ce_loss
            else:
                focal_loss = focal_weight * ce_loss

            # Average loss for this class (ignored pixels already removed by flatten)
            loss_c = focal_loss.mean()
            losses.append(loss_c)

        # Return average loss over both classes
        return torch.stack(losses).mean()

# class TverskyLoss(nn.Module):
    # def __init__(self, ignore_index=255, smooth=1.0, alpha=0.5, beta=0.5):
        # super(TverskyLoss, self).__init__()
        # self.ignore_index = ignore_index
        # self.smooth = smooth
        # self.alpha = alpha
        # self.beta = beta

    # def forward(self, input, target):
        # input, target = flatten(input, target, self.ignore_index)
        # input = softmax(input, dim=1)
        # num_classes = input.size(1)
        # losses = []
        # for c in range(num_classes):
            # target_c = (target == c).float()
            # input_c = input[:, c]
            
            # t_p = (input_c * target_c).sum()
            # f_p = ((1-target_c) * input_c).sum()
            # f_n = (target_c * (1-input_c)).sum()
            # tversky = (t_p + self.smooth) / (t_p + self.alpha*f_p + self.beta*f_n + self.smooth)
            
            # losses.append(1-tversky)
        
        # losses = torch.stack(losses)
        # loss = losses.mean()
        # return loss
        

from scipy.ndimage import distance_transform_edt

class BoundaryLoss(nn.Module):
    """
    基于有符号距离场 (SDF) 的边界损失，专为单通道二分类分割设计。
    输入: pred (B,1,H,W) 经过 Sigmoid 的概率图，值域 [0,1]
         target (B,1,H,W) 二值标签，0/1
    输出: 边界损失标量
    """
    def __init__(self, boundary_width=3):
        super(BoundaryLoss, self).__init__()
        self.boundary_width = boundary_width  # 考虑边界多宽范围内的像素

    def _compute_sdf(self, mask_np):
        """计算单张二值 mask 的有符号距离场 (SDF)"""
        # 前景距离
        pos_dist = distance_transform_edt(mask_np)
        # 背景距离
        neg_dist = distance_transform_edt(1 - mask_np)
        sdf = pos_dist - neg_dist
        # 将 SDF 截断到 [-boundary_width, boundary_width] 并归一化到 [-1,1]
        sdf = np.clip(sdf, -self.boundary_width, self.boundary_width)
        sdf = sdf / self.boundary_width   # 现在范围 [-1,1]
        return sdf

    def forward(self, pred, target):
        """
        pred:   (B,1,H,W) 预测概率，值域 [0,1]
        target: (B,1,H,W) 真实标签，0/1
        """
        batch_size = pred.size(0)
        device = pred.device
        total_loss = 0.0

        # 将 target 转为 numpy 并计算 SDF（这部分不参与梯度，只作为静态权重）
        target_np = target.detach().cpu().numpy()  # (B,1,H,W)
        for i in range(batch_size):
            mask = target_np[i, 0]                 # (H,W)
            sdf = self._compute_sdf(mask)          # (H,W), numpy
            sdf_tensor = torch.from_numpy(sdf).float().to(device)

            # 构造边界权重：离边界越近（|sdf| 越小）权重越大
            # 使用高斯权重，或者简单使用 exp(-sdf^2 / (2*sigma^2))
            sigma = self.boundary_width / 3.0
            weight = torch.exp(- (sdf_tensor ** 2) / (2 * sigma ** 2))

            # 计算加权 BCE 损失
            pred_i = pred[i, 0]      # (H,W)
            target_i = target[i, 0]  # (H,W)
            bce = F.binary_cross_entropy(pred_i, target_i, reduction='none')
            loss_i = (bce * weight).sum() / (weight.sum() + 1e-8)
            total_loss += loss_i

        return total_loss / batch_size


class BoundaryDoULoss(nn.Module):
    """
    Boundary Difference over Union Loss
    完全复现原始实现，但输入为单通道概率图（前景概率），标签为二值 mask。
    内部会构造背景和前景两个类别的 softmax 概率，然后按原始逻辑计算。
    """
    def __init__(self, n_classes=2, smooth=1e-5, alpha_trunc=0.8):
        super(BoundaryDoULoss, self).__init__()
        self.n_classes = n_classes   # 固定为 2
        self.smooth = smooth
        self.alpha_trunc = alpha_trunc
        # 定义边界检测卷积核（十字形）
        kernel = torch.tensor([[0, 1, 0],
                               [1, 1, 1],
                               [0, 1, 0]], dtype=torch.float32).view(1, 1, 3, 3)
        self.register_buffer('kernel', kernel)

    def _one_hot_encoder(self, input_tensor):
        """将 (B,1,H,W) 的二值标签转为 (B,2,H,W) one-hot"""
        tensor_list = []
        for i in range(self.n_classes):
            temp_prob = (input_tensor == i).float()
            tensor_list.append(temp_prob)
        output_tensor = torch.cat(tensor_list, dim=1)
        return output_tensor

    def _adaptive_size(self, score, target):
        """
        score:  (B,H,W) 某一类别的预测概率（连续值，0~1）
        target: (B,H,W) 该类别的 one-hot 二值图
        返回该类的损失分量
        """
        # 边界比例系数 alpha
        # 对 target 做卷积，检测边界
        # 注意 target 是二值 (0/1)，需要加 batch 和 channel 维度
        target_4d = target.unsqueeze(1)  # (B,1,H,W)
        conv = F.conv2d(target_4d, self.kernel.to(target_4d.device), padding=1)  # (B,1,H,W)
        conv = conv.squeeze(1)           # (B,H,W)
        # 只考虑前景像素（target==1）的卷积结果
        conv = conv * target             # 非前景像素置 0
        # 将值为 5 的像素置 0（内部像素，非边界）
        conv[conv == 5] = 0
        C = (conv > 0).sum(dim=(1,2)).float()   # 边界像素数 (B,)
        S = target.sum(dim=(1,2)).float()       # 前景总像素数 (B,)
        alpha = 1 - (C + self.smooth) / (S + self.smooth)
        alpha = 2 * alpha - 1
        alpha = torch.clamp(alpha, max=self.alpha_trunc)  # 截断

        # 计算 DoU 损失
        intersect = (score * target).sum(dim=(1,2))
        y_sum = (target * target).sum(dim=(1,2))
        z_sum = (score * score).sum(dim=(1,2))

        numerator = z_sum + y_sum - 2 * intersect + self.smooth
        denominator = z_sum + y_sum - (1 + alpha) * intersect + self.smooth
        loss = numerator / denominator   # (B,)
        return loss.mean()   # 对 batch 取平均

    def forward(self, pred, target):
        """
        pred:   (B,1,H,W) 预测概率图（前景概率），值域 [0,1]
        target: (B,1,H,W) 二值标签，0 或 1
        """
        # 构造 softmax 概率分布：背景概率 = 1-pred，前景概率 = pred
        # 得到 (B,2,H,W)
        prob_bg = 1 - pred
        prob_fg = pred
        score = torch.cat([prob_bg, prob_fg], dim=1)  # (B,2,H,W)

        # one-hot 编码标签
        target_onehot = self._one_hot_encoder(target)  # (B,2,H,W)

        # 对每个类别计算损失并平均
        total_loss = 0.0
        for c in range(self.n_classes):
            loss_c = self._adaptive_size(score[:, c, :, :], target_onehot[:, c, :, :])
            total_loss += loss_c
        return total_loss / self.n_classes