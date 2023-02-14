from vedo import settings as vsettings
import sys

DEBUG = False  # set to True to see more detailed logs

# ------------------------------- vedo settings ------------------------------ #

vsettings.pointSmoothing = False
vsettings.lineSmoothing = False
vsettings.polygonSmoothing = False
vsettings.immediateRendering = False

vsettings.useDepthPeeling = True
vsettings.alphaBitPlanes = 1
vsettings.maxNumberOfPeels = 12
vsettings.occlusionRatio = 0.1
vsettings.multiSamples = 0 if sys.platform == "darwin" else 8

vsettings.useSSAO = True

# For transparent background with screenshots
vsettings.screenshotTransparentBackground = False  # vedo for transparent bg
vsettings.useFXAA = True  # This needs to be false for transparent bg


# --------------------------- brainrender settings --------------------------- #

BACKGROUND_COLOR = "white"
DEFAULT_ATLAS = "allen_mouse_25um"  # default atlas
DEFAULT_CAMERA = "techpaper_cam_01"  # Default camera settings (orientation etc. see brainrender.camera.py)
INTERACTIVE = True  # rendering interactive ?
LW = 1  # e.g. for silhouettes
ROOT_COLOR = [0.8, 0.8, 0.8]  # color of the overall brain model's actor
ROOT_ALPHA = 0.1  # transparency of the overall brain model's actor'
SCREENSHOT_SCALE = 8
SHADER_STYLE = "plastic"  # affects the look of rendered brain regions: [metallic, plastic, shiny, glossy, cartoon]
SHOW_AXES = False
WHOLE_SCREEN = False  # If true render window is full screen
OFFSCREEN = True

