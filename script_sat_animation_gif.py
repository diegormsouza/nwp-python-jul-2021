#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC Training: NWP Data Processing With Python - Script 10: Creating an Animation (Satellite)
# Author: Diego Souza
#-----------------------------------------------------------------------------------------------------------
import imageio        # Python interface to read and write a wide range of image data

# Images we want to include in the GIF
files = ['Animation/G16_B13_202102181800.png', 'Animation/G16_B13_202102181810.png', 'Animation/G16_B13_202102181820.png', 'Animation/G16_B13_202102181830.png', 'Animation/G16_B13_202102181840.png', 'Animation/G16_B13_202102181850.png', 
         'Animation/G16_B13_202102181900.png', 'Animation/G16_B13_202102181910.png', 'Animation/G16_B13_202102181920.png', 'Animation/G16_B13_202102181930.png', 'Animation/G16_B13_202102181940.png', 'Animation/G16_B13_202102181950.png']

# Create the GIF
images = []
for file in files:
    images.append(imageio.imread(file))
imageio.mimsave('sat_animation.gif', images, fps=1)

# Open the GIF
print("\nOpening the GIF..\n")
from IPython.display import Image
Image(open('sat_animation.gif','rb').read())
