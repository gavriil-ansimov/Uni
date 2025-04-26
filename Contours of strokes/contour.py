#!/usr/bin/env python3
import imutils
import cv2
import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt

image = cv2.imread("D:/Python/13-5.bmp")

image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

gray = image[:, :, 0]

# Noise reduction
#blurred = cv2.GaussianBlur(gray, (5, 5), 0)
#blurred  = cv2.medianBlur(gray, 3)
bblurred = cv2.bilateralFilter(gray, 7, 75, 75)
# thresholding
# thresh = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)[1]
thresh = cv2.threshold(bblurred, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)[1]

se = np.ones((7,7), dtype='uint8')
opening = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, se)
closing = cv2.morphologyEx(opening, cv2.MORPH_CLOSE, se)
#plt.imshow(closing,'gray')
#plt.show()
edge = cv2.Canny(closing, 100, 200)

contours, h = cv2.findContours(image=edge, mode=cv2.RETR_CCOMP, method=cv2.CHAIN_APPROX_SIMPLE)
contours = sorted(contours, key=cv2.contourArea, reverse=True) 
image_copy = image.copy()
#cnts = imutils.grab_contours(contours)

k = 0
area = -1
x1 = -100
c = [x1]
for i in contours:
    flag = True
    M = cv2.moments(i)
    # Area of contour
    area = M["m00"]
    # Centroid of contour
    x = int(M["m10"] / M["m00"])
    y = int(M["m01"] / M["m00"])
    #print(c)
    for j in range(len(c)):
        if (x > c[j] - 30 and x < c[j] + 30):
            flag = False
    image_copy = image.copy()
    lx = len(image_copy)
    ly = len(image_copy[1])
    sum1, sum2, sum3 = 0, 0, 0
    cnt = 0 
    if flag:
        cv2.drawContours(image=image_copy, contours=[i], contourIdx=-1, color=(0, 255, 0), thickness=cv2.FILLED, lineType=cv2.LINE_AA)
        for a in range(lx):
            for b in range(ly):
                if image_copy[a][b][1] == 255:
                    sum1 += image[a][b][0]
                    sum2 += image[a][b][1]
                    sum3 += image[a][b][2]
                    cnt += 1
        print(sum1/cnt, sum2/cnt, sum3/cnt)
        cv2.imwrite(f'res{k}.jpg', image_copy)
        k += 1
        c.append(x)
    if k >= 5:
        break
