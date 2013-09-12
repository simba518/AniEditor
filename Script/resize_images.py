import glob
import os
import Image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# more informations of image processing using python could be found at:
# http://www.riisen.dk/dop/pil.html
# http://www.pythonware.com/library/pil/handbook/image.htm
# http://matplotlib.sourceforge.net/users/image_tutorial.html

# resize all the images in the "src_image_dir" directory with the size specified
# by the "box", then save it to the target_image_dir with the old_name replaced
# by new_name.
def resize_images_cname (src_image_dir, target_image_dir, box, old_name, new_name , appendix = "png"):
    # get the filenames of the src directory and the target directory.
    src_image_files = []
    target_image_files = []
    os.chdir(src_image_dir)
    for files in glob.glob("*."+appendix):
        src_image_files.append( os.path.join(src_image_dir,files) )
        target_image_files.append( os.path.join(target_image_dir,files) )

    # read in all files in a directory
    images = []
    for f in src_image_files:
        images.append(Image.open(f) )

    # resize all images
    new_images = []
    for img in images:
        new_images.append( img.crop(box) )
        
    # save all images to the target directory
    if not os.path.exists(target_image_dir):
        os.mkdir(target_image_dir)

    for i in range(0,len(target_image_files)):
        new_images[i].save(target_image_files[i].replace(old_name,new_name))



# resize all the images in the "src_image_dir" directory with the size specified
# by the "box", then save it to the target_image_dir with the same name.
def resize_images (src_image_dir, target_image_dir, box, appendix = "png"):
    resize_images_cname(src_image_dir, target_image_dir, box, "","",appendix)

# show an image
def show_image(image_name):
    img = mpimg.imread(image_name)
    imgplot = plt.imshow(img)
    plt.show()
