#----------------------------------------------------------------------------------------------------------------------
# INPE / CPTEC Training: NWP Data Processing With Python - Script 10: Creating an Animation
# Author: Diego Souza
#----------------------------------------------------------------------------------------------------------------------
import imageio        # Python interface to read and write a wide range of image data

# Images we want to include in the GIF
files = ['image_loop_0.png', 'image_loop_3.png', 'image_loop_6.png', 'image_loop_9.png', 'image_loop_12.png', 
         'image_loop_15.png', 'image_loop_18.png', 'image_loop_21.png', 'image_loop_24.png']

# Create the GIF
images = []
for file in files:
    images.append(imageio.imread('Animation//' + file))

# Save the GIF
imageio.mimsave('Animation//animation.gif', images, fps=1)