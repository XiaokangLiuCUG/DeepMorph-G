import os
import pandas as pd
from skimage import io, transform
import torch
import torchvision
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms#, utils
# import torch.optim as optim
import csv
import numpy as np
from PIL import Image
import glob

from data_loader import RescaleT
from data_loader import ToTensor
from data_loader import ToTensorLab
from data_loader import SalObjDataset

from model import U2NET # full size version 173.6 MB
from model import U2NETP # small version u2net 4.7 MB

def get_image_paths_from_excel(excel_path, image_root, sheet_name='Sheet1'):
    if not os.path.exists(excel_path):
        raise FileNotFoundError(f"Excel file does not exist: {excel_path}")
    
    df = pd.read_excel(excel_path, sheet_name=sheet_name)
    
    if 'file_name' not in df.columns:
        raise KeyError(f"'file_name' column not found in sheet '{sheet_name}'. Existing columns: {df.columns.tolist()}")
    
    # Drop missing values and empty strings, and convert to string type
    file_names = df['file_name'].dropna().astype(str).str.strip()
    file_names = file_names[file_names != '']
    
    abs_paths = []
    missing_count = 0
    for fname in file_names:
        full_path = os.path.join(image_root, fname)
        if os.path.exists(full_path):
            abs_paths.append(full_path)
        else:
            missing_count += 1
            # Optional: print a warning for missing files (uncomment to see specific missing files)
            # print(f"Warning: image does not exist - {full_path}")
    
    if missing_count > 0:
        print(f"A total of {missing_count} image files do not exist and were skipped.")
    print(f"Number of valid images: {len(abs_paths)} / {len(file_names)}")
    
    return abs_paths


                                       
# normalize the predicted SOD probability map
def normPRED(d):
    ma = torch.max(d)
    mi = torch.min(d)

    dn = (d-mi)/(ma-mi)

    return dn

def save_output(image_name,pred,d_dir):

    predict = pred
    predict = predict.squeeze()
    predict_np = predict.cpu().data.numpy()

    im = Image.fromarray(predict_np*255).convert('RGB')
    img_name = image_name.split(os.sep)[-1]
    image = io.imread(image_name)
    imo = im.resize((image.shape[1],image.shape[0]),resample=Image.BILINEAR)
    #convert to binary images
    pb_np = np.array(imo)
    thresh = 200
    fn = lambda x: 255 if x < thresh else 0
    r = imo.point(fn)#, mode='1')
    
    aaa = img_name.split(".")
    bbb = aaa[0:-1]
    imidx = bbb[0]
    for i in range(1,len(bbb)):
        imidx = imidx + "." + bbb[i]

    #imo.save(d_dir+imidx+'.png')
    r.save(d_dir+imidx+'.png')
    
def main():

    # --------- 1. get image path and name ---------
    model_name='u2net'#u2netp
    file_name = model_name

    image_dir = r'D:\项目\DeepMorph-G\s fig1' 
    image_extensions = ('.jpg', '.jpeg', '.png', '.bmp', '.tif', '.tiff')

    img_name_list = []
    for root, dirs, files in os.walk(image_dir):
        for f in files:
            if f.lower().endswith(image_extensions):
                img_name_list.append(os.path.join(root, f))

    print(img_name_list[0])
    print(len(img_name_list))
    clade = 'three'#aperture outline  GBIFgastropoda
    prediction_dir = os.path.join(os.getcwd(), 'test_data', clade + os.sep)
    
    #model_dir = os.path.join(os.getcwd(), 'saved_models', model_name, file_name + '.pth')
    model_dir = r'E:\Python_learning\PyTorch\U-2-Net-master\saved_models_fossil\u2net\u2net_outline_jaccardloss+0.5=loss_boundary_max_lr=0.005_div_factor=5.0, final_div_factor=10_best_class1_iou.pth'#E:\Python_learning\PyTorch\U-2-Net-master\saved_models\u2net\u2net.pth'

    # --------- 2. dataloader ---------
    #1. dataloader
    test_salobj_dataset = SalObjDataset(img_name_list = img_name_list,
                                        lbl_name_list = [],
                                        transform=transforms.Compose([RescaleT(320),#Highet with
                                                                      ToTensorLab(flag=0)])
                                        )
    test_salobj_dataloader = DataLoader(test_salobj_dataset,
                              
                              batch_size=1,
                                        shuffle=False,
                                        num_workers=0)

    # --------- 3. model define ---------
    if(model_name=='u2net'):
        print("...load U2NET---173.6 MB")
        net = U2NET(3,1)
    elif(model_name=='u2netp'):
        print("...load U2NEP---4.7 MB")
        net = U2NETP(3,1)

    if torch.cuda.is_available():
        net.load_state_dict(torch.load(model_dir))
        net.cuda()
    else:
        net.load_state_dict(torch.load(model_dir, map_location='cpu'))
    net.eval()

    ##########################print model
    from torchsummary import summary
    summary(net, (3, 320, 320))
    ##########################
    # --------- 4. inference for each image ---------
    for i_test, data_test in enumerate(test_salobj_dataloader):

        print("inferencing:",img_name_list[i_test].split(os.sep)[-1])

        inputs_test = data_test['image']
        inputs_test = inputs_test.type(torch.FloatTensor)

        if torch.cuda.is_available():
            inputs_test = Variable(inputs_test.cuda())
        else:
            inputs_test = Variable(inputs_test)

        d1,d2,d3,d4,d5,d6,d7= net(inputs_test)

        # normalization
        pred = d1[:,0,:,:]
        pred = normPRED(pred)

        # save results to test_results folder
        if not os.path.exists(prediction_dir):
            os.makedirs(prediction_dir, exist_ok=True)
        save_output(img_name_list[i_test],pred,prediction_dir)

        del d1,d2,d3,d4,d5,d6,d7

if __name__ == "__main__":
    main()