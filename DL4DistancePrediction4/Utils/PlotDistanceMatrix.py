import os
import shutil
import commands
import sys
import numpy as np
from PIL import Image, ImageDraw

## all distance smaller than the below cutoff will be displayed as black
MinDist4Display=3.0

# all distance larger than the below cutoff will be displayed as white
MaxDist4Display=15.0

NUM_HEIGHT = 12
FINAL_SIZE = 624
NUMS_TO_SHOW = 9.5625

#------- basic command for processing figures --------# start
def GetColor(cont):
    #return (255*(1-cont), 255*(1-cont), 255*(1-cont))
    return (int(255*(1-cont)), int(255*(1-cont)), int(255*(1-cont)))

def FillContacts(contacts, cont, i, j, multiplier):
    for k in range(0,multiplier):
        for l in range(0,multiplier):
            contacts[(i-1)*multiplier + k][(j-1)*multiplier + l] = cont
            contacts[(j-1)*multiplier + k][(i-1)*multiplier + l] = cont


def GetPixels(cfile, multiplier, lenseq):
    ln_new = ''
    ln_old = ''
    while 1:
        ln_old = ln_new
        ln_new = cfile.readline()
        parts = ln_new.split()
        if len(parts) > 0:
            if parts[0].isdigit():
                break

    contacts = []

    for i in range(0, lenseq*multiplier):
        contacts.append([])
        for j in range(0, lenseq*multiplier):
            contacts[i].append(0.0)
    if len(parts) >= 3:
        i = int(parts[0].strip())
        j = int(parts[1].strip())
        cont = float(parts[-1].strip())
        FillContacts(contacts, cont, i, j, multiplier)

    for ln in cfile:
        parts = ln.split()
        if len(parts) >= 3:
            i = int(parts[0].strip())
            j = int(parts[1].strip())
            cont = float(parts[-1].strip())
            FillContacts(contacts, cont, i, j, multiplier)

    #bordersize = int((NUM_HEIGHT*FINAL_SIZE)/(lenseq*multiplier))+1
    bordersize = int((NUM_HEIGHT+FINAL_SIZE)/(lenseq*multiplier))
    if bordersize < 1:
        bordersize = 1
    colors = []
    for i in range(0, lenseq*multiplier + bordersize):
        for j in range(0, bordersize):
            colors.append((255,255,255))
    for i in contacts:
        for j in range(0, bordersize):
            colors.append((255,255,255))
        for j in i:
            colors.append(GetColor(j))

    return (lenseq, colors)


def CreateSequenceNums(imglist, imglens):
    img_hor = Image.new('RGB', (FINAL_SIZE-NUM_HEIGHT, NUM_HEIGHT))
    img_temp = Image.new('RGB', (FINAL_SIZE-NUM_HEIGHT, NUM_HEIGHT))
    colors = [(0,0,0)]*(FINAL_SIZE-NUM_HEIGHT)*NUM_HEIGHT

    for i in range(9*(FINAL_SIZE-NUM_HEIGHT), 10*(FINAL_SIZE-NUM_HEIGHT)):
        colors[i] = (255,255,255)

    img_hor.putdata(colors)
    img_temp.putdata(colors)
    newspacing = (FINAL_SIZE-NUM_HEIGHT)/NUMS_TO_SHOW

    for i in range(0, 10):
        #img_hor.paste(imglist[i], (round(i*newspacing),0))
        #img_temp.paste(imglist[i], (FINAL_SIZE-NUM_HEIGHT-round(i*newspacing)-imglens[i],0))

	img_hor.paste(imglist[i], (int(round(i*newspacing)),0))
        img_temp.paste(imglist[i], (FINAL_SIZE-NUM_HEIGHT-int(round(i*newspacing))-imglens[i],0))

    img_ver = img_temp.transpose(Image.ROTATE_90)

    return (img_hor, img_ver)

#------- basic command for processing figures --------# over



#---------- main --------#
if len(sys.argv) < 2:
	print 'python PlotDistanceMatrix.py targetName.bound.txt [result_dir]'
	print '\ttargetName.bound.txt: a file for the predicted contact distance matrix in text format'
	print '\tresult_dir: optional, saving the resultant images in this folder; default is current work directory'
    	exit(-1)

contactFile = sys.argv[1]
result_dir="./"
targetName=os.path.basename(contactFile).split('.')[0]

if len(sys.argv)>=3:
	result_dir = sys.argv[2]


if not os.path.isfile(contactFile):
	print 'the contact distance file (in txt format) does not exist: ', contactFile
	exit(-1)



lenseq=0   ## seq len
MatrixDistance=[]  ## coordinate data


#---------- read file, get data --------#
with open(contactFile, 'r') as f:
    for one in f.readlines():
        if one.startswith("SEQ"):          ## the seq
            list=one.strip().split(" ")
            lenseq=lenseq+len(list[1])
        if one[0].isdigit():               ## read data
            list=one.split(' ')
            MatrixDistance.append([list[0],list[1],list[2]])
        if one.startswith("END"):  break  ## the end
#print(lenseq,MatrixDistance)




#---------- analyze data, get matrix --------#



conts = []

for i in xrange(lenseq):
	onerow = [0.]*lenseq
	conts.append(onerow)

for ln in MatrixDistance:
	i = int(ln[0])
	j = int(ln[1])
	if i > lenseq or j > lenseq:
		print 'residue indices in this row are out of range: ', ln
		exit(-1)

        cont = np.float32(ln[2])
        if cont <= MinDist4Display : cont = 1    ##  black
        elif cont >= MaxDist4Display : cont= 0    ##  white
        else : cont= (MaxDist4Display - cont) / (MaxDist4Display - MinDist4Display) 
        conts[i - 1][j - 1] = cont
        conts[j - 1][i - 1] = cont



## convert probability values to RGB, note that colors is similar to just a 1D array.
colors = []
for x in conts:
	for y in x:
        	colors.append((int(255 * (1 - y)), int(255 * (1 - y)), int(255 * (1 - y))))

img1 = Image.new('RGB', (lenseq, lenseq))
img1.putdata(colors)
img2 = img1.resize((FINAL_SIZE - NUM_HEIGHT, FINAL_SIZE - NUM_HEIGHT))

#img2.save(os.path.join(result_dir,  targetName + 'img2.png') )

## show residue numbers
mynum = lenseq / NUMS_TO_SHOW  # at which interval to space numbers if we are displaying 10(each digit is 12 pixels)
imglist = []
imglens = []
for i in range(0, 10):
    	indx = int(i * mynum) + 1
    	mytext = str(indx)
    	if indx < 10:
        	imglen = 1 * NUM_HEIGHT / 2
    	elif indx < 100:
        	imglen = 2 * NUM_HEIGHT / 2
    	else:
        	imglen = 3 * NUM_HEIGHT / 2
    	imglens.append(imglen)
    	imgtemp = Image.new('RGB', (imglen, NUM_HEIGHT))
    	imgDrawer = ImageDraw.Draw(imgtemp)
    	imgDrawer.text((0, 0), mytext)
    	imglist.append(imgtemp)

imgnums = CreateSequenceNums(imglist, imglens)
img3 = Image.new('RGB', (FINAL_SIZE, FINAL_SIZE))
img3.paste(imgnums[0], (NUM_HEIGHT, 0))
img3.paste(imgnums[1], (0, NUM_HEIGHT))
img3.paste(img2, (NUM_HEIGHT, NUM_HEIGHT))
img4 = Image.new('RGB', (NUM_HEIGHT, NUM_HEIGHT))
coWhite = [(255, 255, 255)] * NUM_HEIGHT * NUM_HEIGHT
img4.putdata(coWhite)
img3.paste(img4, (0, 0))
img3 = img3.resize((FINAL_SIZE, FINAL_SIZE))
img3.save(os.path.join(result_dir,  targetName + '_dist.png') )

## generate multiple images of different size
for num in range(1,6):
    	img4=img3.resize((FINAL_SIZE*num, FINAL_SIZE*num))
	img4.save(os.path.join(result_dir,  targetName + '_' + str(num)+ '_dist.png') )
