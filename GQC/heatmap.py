import sys
import importlib.util
import pybedtools
import logging

if (spec := importlib.util.find_spec('matplotlib')) is not None:
    #importing the present module
    module = importlib.util.module_from_spec(spec)
    sys.modules['matplotlib'] = module
    spec.loader.exec_module(module)
    matlibimported = True
else:
    print(f"To print heatmap colors in bedfiles, you need to install matplotlib")
    matlibimported = False

if matlibimported:
   import matplotlib.pyplot as plt
   import matplotlib.colors as mcolors

logger = logging.getLogger(__name__)

def add_colors_to_intervals(windowedintervalswithcounts, args):
    if not matlibimported:
        return None

    intervalbedstringwithcolors = ""

    allscores = []
    for interval in windowedintervalswithcounts:
        score = int(interval[4])
        if args.windowsize > interval.end - interval.start:
            score = int(score * args.windowsize / (interval.end - interval.start))
        allscores.append(score)

    if args.heatmapmincount and args.heatmapmaxcount:
        minscore = args.heatmapmincount
        maxscore = args.heatmapmaxcount
    else:
        numintervals = len(allscores)
        lowcutoff = int(0.05*numintervals)
        highcutoff = int(0.95*numintervals)
        allscores.sort()
        minscore = allscores[lowcutoff]
        maxscore = allscores[highcutoff]
    logger.info("Using count " + str(minscore) + " for minimum end of spectrum and count " + str(maxscore) + " for maximum end of spectrum")

    for interval in windowedintervalswithcounts:
        ref = interval.chrom
        start = interval.start
        end = interval.end
        name = interval.name
        score = int(interval[4])

        rgbval = convert_score_to_rgb(score, min_score=minscore, max_score=maxscore)

        intervalbedstringwithcolors = intervalbedstringwithcolors + ref + "\t" + str(start) + "\t" + str(end) + "\t" + name + "\t" + str(score) + "\t.\t" + str(start) + "\t" + str(end) + "\t" + rgbval + "\n"

    intervalbedtool =  pybedtools.BedTool(intervalbedstringwithcolors, from_string=True)

    return intervalbedtool

# Function to convert a score to an RGB color (adapted from a function written by Steven Solar)
def convert_score_to_rgb(score:int, min_score=0, max_score=1000, colormap='YlOrRd')->str:
    if score < min_score:
        score = min_score
    elif score > max_score:
        score = max_score

    # Define the colormap
    cmap = plt.get_cmap(colormap)

    # Define the normalization range
    norm = mcolors.Normalize(vmin=min_score, vmax=max_score)

    # Get the RGBA color corresponding to the percent identity
    rgba = cmap(norm(score))

    # Convert RGBA to RGB
    rgb = tuple(int(255 * x) for x in rgba[:3])
    return str(rgb[0]) + "," + str(rgb[1]) + "," + str(rgb[2])


