---
layout: tutorial_hands_on

title: Quantitative Analysis of Histological Staining Using Color Deconvolution
zenodo_link: ''
questions:
- How can I quantify the percentage of stained area in histological images?
- How does color deconvolution separate individual stain components from brightfield images?
- How can I apply this workflow to IHC and Masson's Trichrome stained tissue sections?
objectives:
- Apply color deconvolution to separate stain channels in histological images
- Extract and isolate the stain channel of interest (e.g. DAB, aniline blue)
- Apply automatic thresholding to distinguish stained from unstained regions
- Calculate the percentage of positively stained area relative to total tissue area
- Interpret quantitative staining results across multiple images
time_estimation: 3H
key_points:
- Color deconvolution mathematically separates RGB images into individual stain components based on known absorption spectra
- The same H-E-DAB (HED) preset is used for both IHC and Masson's Trichrome; channel selection determines the stain of interest
- Otsu thresholding provides an automated, reproducible way to distinguish stained from unstained pixels
- The percent stained area is calculated as (stained pixels / total tissue pixels) × 100
- This workflow accepts a collection of images for batch quantification across multiple samples
contributors:
- dianichj
---

Manually scoring histological staining across dozens of images is time-consuming and subjective. Two researchers looking at the same slide may reach different conclusions about how much staining is present. Computational quantification solves this: it applies the same criteria to every image, produces a numeric result, and scales to large datasets without additional effort.

This tutorial walks you through a Galaxy workflow that quantifies stained area in brightfield histological images — from raw microscopy image to a final table of percentages ready for statistical analysis. The approach is built around **color deconvolution**, a technique that mathematically separates overlapping stain signals so you can measure each one independently.

You will work with two types of histological staining:

- **IHC (Immunohistochemistry)** — detecting CD11b-positive myeloid cells using a DAB chromogen
- **Masson's Trichrome (MT)** — quantifying collagen deposition as a measure of cardiac fibrosis

Both use the same core workflow; the only difference is which channel you extract after deconvolution.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background: What Is Color Deconvolution and Why Do We Need It?

When two stains are applied to the same tissue section, for example, hematoxylin (blue) and DAB (brown) in IHC, or Weigert's hematoxylin, Biebrich scarlet, and aniline blue in Masson's Trichrome, their colors overlap in the RGB image. This means that one cannot simply threshold on "brownness" or "blueness" because the RGB channels mix all stain signals together.

**Color deconvolution** solves this by using a **stain matrix**: a set of vectors that describe how much each stain absorbs light in the red, green, and blue channels. The algorithm inverts this matrix to compute the optical density of each stain independently at every pixel, producing one grayscale image per stain component. Higher pixel values in the output correspond to stronger staining.

In this tutorial, we use the **H-E-DAB (HED)** preset for both staining types:

| Channel | IHC interpretation | Masson's Trichrome interpretation |
|---|---|---|
| Channel 1 (index 0) | Hematoxylin | Weigert's hematoxylin (nuclei) |
| Channel 2 (index 1) | **DAB — stain of interest** | Biebrich scarlet (muscle/cytoplasm) |
| Channel 3 (index 2) | Residual | **Aniline blue — stain of interest** |

> <details-title>What is a stain vector?</details-title>
>
> A stain vector describes how much a given stain absorbs light in each of the three RGB color channels. For example, DAB absorbs strongly in the blue channel and weakly in the red channel, while hematoxylin absorbs more in the red and green channels. These values are determined empirically from pure stain reference images or taken from published standards. The color deconvolution algorithm uses these vectors to solve a system of linear equations and separate the mixed stain signals pixel by pixel.
>
{: .details}

> <details-title>When to consider Non-negative Matrix Factorization (NMF) instead</details-title>
>
> For most standard IHC images with a clean DAB + hematoxylin combination, the HED preset is sufficient. However, when staining is uneven, the signal is weak, or there is significant spectral overlap, a data-driven approach may give better results.
>
> **Non-negative Matrix Factorization (NMF)** learns the stain components directly from the image data, without relying on predefined stain vectors. This makes it more flexible when staining deviates from standard reference spectra — for example, due to differences in staining protocols, tissue processing, or scanner calibration. To use NMF, select `Non-negative matrix factorization` as the transformation type. Since NMF components are not labeled by stain name, you will need to visually inspect the output channels to identify which one corresponds to your stain of interest.
>
{: .details}

# The Dataset: Cardiac Tissue After Myocardial Infarction

The images in this tutorial come from a study investigating the role of 4-oxo retinoic acid (4-oxo RA) in maintaining hematopoietic stem cell dormancy after myocardial infarction (MI) {% cite Rettkowski2025 %}. Understanding this biological context helps interpret what the numbers mean.

**The experimental model:** Mouse hearts were subjected to LAD coronary artery ligation to induce MI. Animals were treated with 4-oxo RA or with the vehicle solution (DSMO) after MI. Serial cardiac sections were stained to assess two aspects of the post-infarction response:

**IHC for CD11b** detects myeloid leukocytes (monocytes, macrophages, neutrophils) — a readout of local inflammation. Higher CD11b-positive area = more immune cell infiltration after MI.

**Masson's Trichrome** stains collagen blue, allowing quantification of fibrotic remodeling. Higher collagen area = more fibrosis, reflecting more extensive tissue damage.

The hypothesis is that 4-oxo RA reduces immune cell mobilization from bone marrow, leading to less inflammatory infiltration (reduced presence of CD11b+ cells), less fibrosis, and better cardiac tissue preservation post-MI.

**What you are working with:** In the original study, regions of interest (ROIs) were manually selected from whole-slide images to focus on specific anatomical zones (Infarct, Border, Remote) and exclude tissue artifacts. For this tutorial, those ROIs have already been cropped and are provided as ready-to-use image collections — so you can focus entirely on the quantification workflow.

Here is an example of what the raw input images look like:

![Example IHC image showing CD11b staining in brown (DAB) with hematoxylin counterstain in blue. The brown signal marks myeloid cell infiltration in cardiac tissue post-MI.](../../images/example_IHC_input.png "Example input: IHC for CD11b in cardiac tissue after myocardial infarction.")

![Example Masson's Trichrome image showing collagen in blue (aniline blue), muscle in red (Biebrich scarlet), and nuclei in dark purple. Blue area reflects fibrotic remodeling.](../../images/example_MT_input.png "Example input: Masson's Trichrome staining in cardiac tissue after myocardial infarction.")

> <comment-title>Sample data vs. full dataset</comment-title>
>
> The images provided here are a representative subset of the full dataset from {% cite Rettkowski2025 %}. In practice, workflows like this one are run across large image batches spanning multiple experiments, animals, conditions, and anatomical zones. The workflow you will run here is similar to what was applied at scale in the original study.
>
{: .comment}

## Data Upload

> <hands-on-title>Upload your images</hands-on-title>
>
> 1. Create a new history and name it something meaningful, e.g. `Histology Stain Quantification`
> 2. Import the images from [Zenodo]({{ page.zenodo_link }}) or from the shared data library (`GTN - Material` → `{{ page.topic_name }}` → `{{ page.title }}`):
>
>    ```
>    
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets with descriptive names (e.g. `IHC_sample1.tiff`, `MT_sample1.tiff`)
> 4. Check that the datatype is `tiff`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tiff" %}
>
> 5. Organize your images into a **dataset collection** (one collection per staining type). Collections allow the workflow to process all images in a single run.
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
> > <comment-title>Image requirements</comment-title>
> >
> > Images must be brightfield RGB microscopy images in TIFF format. Fluorescence images are not compatible with color deconvolution. Avoid images with strong artifacts, out-of-focus regions, or uneven illumination, as these reduce quantification accuracy.
> >
> {: .comment}
>
{: .hands_on}

{% include _includes/cyoa-choices.html option1="IHC" option2="MT" default="IHC"
   text="This tutorial can be followed with two different staining types: IHC for CD11b, or Masson's Trichrome for collagen quantification. Choose the one that matches your data — the steps are the same, only the channel you extract differs." %}

# Step 1 — Color Deconvolution

The first step separates the mixed stain signals in your brightfield image into individual channels. We use the **H-E-DAB (HED)** color space for both IHC and Masson's Trichrome.

> <hands-on-title>Run color deconvolution</hands-on-title>
>
> 1. {% tool [Perform color deconvolution or transformation](toolshed.g2.bx.psu.edu/repos/imgteam/color_deconvolution/ip_color_deconvolution/0.9+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input image"*: your image collection
>    - *"Transformation type"*: `Deconvolve RGB into Hematoxylin + Eosin + DAB`
>
>    > <comment-title>Output</comment-title>
>    >
>    > Each input image produces one multi-channel TIFF. Each channel is a grayscale image where brighter pixels indicate higher optical density (stronger staining) for that component. You will extract the relevant channel in the next step.
>    >
>    {: .comment}
>
{: .hands_on}

<div class="IHC" markdown="1">

The figure below shows what to expect after deconvolution for an IHC image. The DAB channel (right) clearly isolates the brown signal from the blue hematoxylin counterstain:

![Color deconvolution output for IHC. Left: original RGB image. Middle: Hematoxylin channel. Right: DAB channel isolating the brown CD11b signal.](../../images/deconvolution_output_IHC.png "Color deconvolution output for IHC (CD11b/DAB).")

</div>

<div class="MT" markdown="1">

The figure below shows what to expect after deconvolution for a Masson's Trichrome image. Channel 3 (right) isolates the aniline blue collagen signal:

![Color deconvolution output for Masson's Trichrome. Left: original RGB image. Channels 1–3 shown separately, with channel 3 capturing the aniline blue collagen signal.](../../images/deconvolution_output_MT.png "Color deconvolution output for Masson's Trichrome.")

</div>

# Step 2 — Split Channels and Extract the Stain of Interest

The deconvolution output is a multi-channel TIFF. You need to split it into individual channel images and then pick the one that corresponds to your stain.

> <hands-on-title>Split the multi-channel image</hands-on-title>
>
> 1. {% tool [Split image along axes](toolshed.g2.bx.psu.edu/repos/imgteam/split_image/ip_split_image/2.3.5+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Image to split"*: output of **Perform color deconvolution or transformation** {% icon tool %}
>
{: .hands_on}

This produces a collection of single-channel grayscale images — one per stain component. Now extract the channel you need:

<div class="IHC" markdown="1">

> <hands-on-title>Extract the DAB channel (IHC)</hands-on-title>
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-collection %} *"Input List"*: output of **Split image along axes** {% icon tool %}
>    - *"How should a dataset be selected?"*: `Select by index`
>    - *"Element index"*: `1`
>
>    > <comment-title>Why index 1?</comment-title>
>    >
>    > The split collection is zero-indexed. In the HED preset, Channel 2 (DAB) corresponds to index `1`.
>    >
>    {: .comment}
>
{: .hands_on}

</div>

<div class="MT" markdown="1">

> <hands-on-title>Extract the collagen channel (MT)</hands-on-title>
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-collection %} *"Input List"*: output of **Split image along axes** {% icon tool %}
>    - *"How should a dataset be selected?"*: `Select by index`
>    - *"Element index"*: `2`
>
>    > <comment-title>Why index 2?</comment-title>
>    >
>    > The split collection is zero-indexed. In the HED preset, Channel 3 (aniline blue / collagen) corresponds to index `2`.
>    >
>    {: .comment}
>
{: .hands_on}

</div>

> <question-title></question-title>
>
> 1. What does a high pixel intensity value mean in the deconvolved DAB channel?
> 2. What index would you use to extract the hematoxylin channel from an HED deconvolution?
>
> > <solution-title></solution-title>
> >
> > 1. Higher pixel values in the deconvolved channel represent stronger staining (higher optical density of that dye at that location). Darker brown regions in the original IHC image become bright pixels in the DAB channel.
> > 2. Hematoxylin is Channel 1, so you would use index `0`.
> >
> {: .solution}
>
{: .question}

# Step 3 — Capture Total Image Area

To calculate the *percentage* of stained area, you need two numbers: how many pixels are stained, and how many pixels make up the whole tissue image. You will measure the stained pixels later from the thresholded image. Here you extract the image dimensions to compute the total area.

> <hands-on-title>Get image dimensions from the original image</hands-on-title>
>
> 1. {% tool [Show image info](toolshed.g2.bx.psu.edu/repos/imgteam/image_info/ip_imageinfo/5.7.1+galaxy1) %} with the following parameters:
>    - {% icon param-collection %} *"Input Image"*: your original image collection (the same collection you provided as input to color deconvolution)
>
>    > <comment-title>Why use the original image here?</comment-title>
>    >
>    > We extract dimensions from the raw input image — not the deconvolved output — because the original image faithfully represents the full tissue area captured by the microscope. This ensures the total area denominator is correct regardless of any processing applied downstream.
>    >
>    {: .comment}
>
> 2. {% tool [Select](Grep1) %} to extract the image width:
>    - {% icon param-collection %} *"Select lines from"*: output of **Show image info** {% icon tool %}
>    - *"the pattern"*: `Width =`
>
> 3. {% tool [Select](Grep1) %} to extract the image height:
>    - {% icon param-collection %} *"Select lines from"*: output of **Show image info** {% icon tool %}
>    - *"the pattern"*: `Height =`
>
> 4. {% tool [Text transformation](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sed_tool/9.5+galaxy3) %} to isolate the width value:
>    - {% icon param-collection %} *"File to process"*: output of **Select** (Width) {% icon tool %}
>    - *"SED Program"*: `s/.*= //`
>
>    > <comment-title>What does this expression do?</comment-title>
>    >
>    > The image info tool returns lines like `Width = 1024`. The sed expression `s/.*= //` strips everything up to and including `= `, leaving just the number. This is necessary so it can be used in arithmetic downstream.
>    >
>    {: .comment}
>
> 5. {% tool [Text transformation](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sed_tool/9.5+galaxy3) %} to isolate the height value:
>    - {% icon param-collection %} *"File to process"*: output of **Select** (Height) {% icon tool %}
>    - *"SED Program"*: `s/.*= //`
>
> 6. {% tool [Paste](Paste1) %} to combine width and height side by side:
>    - {% icon param-collection %} *"Paste"*: output of **Text transformation** (width) {% icon tool %}
>    - {% icon param-collection %} *"and"*: output of **Text transformation** (height) {% icon tool %}
>
> 7. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.1+galaxy0) %} to calculate total pixel area (width × height):
>    - {% icon param-collection %} *"Input file"*: output of **Paste** {% icon tool %}
>    - *"Input has a header line with column names?"*: `No`
>    - In *"Expressions"* → *"Insert Expressions"*:
>      - *"Add expression"*: `c1 * c2`
>    - *"If an expression cannot be computed for a row"*: `Fail the entire tool run`
>
{: .hands_on}

> <question-title></question-title>
>
> An image is 1024 × 768 pixels. What is the total pixel area, and why does it matter?
>
> > <solution-title></solution-title>
> >
> > Total pixel area = 1024 × 768 = 786,432 pixels. This is the denominator in the percentage formula: (stained pixels / total pixels) × 100. Without it you can count stained pixels in absolute terms but cannot compare meaningfully across images of different sizes.
> >
> {: .solution}
>
{: .question}

# Step 4 — Threshold the Stain Channel

Now you will convert the extracted grayscale channel into a binary mask: pixels are classified as either stained (value = 1) or unstained (value = 0). This is the foundation for measuring stained area.

> <hands-on-title>Apply Otsu thresholding</hands-on-title>
>
> 1. {% tool [Threshold image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_auto_threshold/ip_threshold/0.25.2+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input image"*: output of **Extract dataset** {% icon tool %} (your extracted stain channel)
>    - *"Thresholding method"*: `Globally adaptive / Otsu`
>
{: .hands_on}

Otsu's method automatically finds the threshold that best separates the two pixel populations — stained and unstained — by minimizing the variance within each group. Because it adapts to each image's intensity distribution, you do not need to set a manual value per image, making the results more consistent across a large batch.

<div class="IHC" markdown="1">

![Thresholding output for IHC. Left: deconvolved DAB channel. Right: binary Otsu mask where white pixels represent DAB-positive area.](../../images/threshold_output_IHC.png "Otsu thresholding output for IHC (CD11b/DAB).")

</div>

<div class="MT" markdown="1">

![Thresholding output for Masson's Trichrome. Left: deconvolved collagen channel. Right: binary Otsu mask where white pixels represent collagen-positive area.](../../images/threshold_output_MT.png "Otsu thresholding output for Masson's Trichrome.")

</div>

> <question-title></question-title>
>
> When might Otsu thresholding give unreliable results?
>
> > <solution-title></solution-title>
> >
> > Otsu works best when the pixel intensity histogram has two distinct peaks (bimodal distribution). It may underperform when staining is very sparse or very dense (producing a unimodal histogram), or when there is significant background noise or uneven illumination. In those cases, consider preprocessing steps such as background subtraction, or apply a manually set threshold offset.
> >
> {: .solution}
>
{: .question}

# Step 5 — Generate ROIs for Visual Validation

Before moving to quantification, it is good practice to visually verify that the thresholded mask captures what you expect. This step detects the stained blobs in the binary mask and generates polygon ROIs around them, which can be inspected or uploaded to OMERO for review alongside the original images.

> <hands-on-title>Detect stained regions and generate ROIs</hands-on-title>
>
> 1. {% tool [Analyze particles](toolshed.g2.bx.psu.edu/repos/imgteam/imagej2_analyze_particles_binary/imagej2_analyze_particles_binary/20240614+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Binary image"*: output of **Threshold image** {% icon tool %}
>    - *"Minimum size"*: `50`
>    - *"Black background"*: `Yes`
>    - *"Export ROIs"*: `Yes`
>
>    > <comment-title>What does the minimum size filter do?</comment-title>
>    >
>    > Setting a minimum particle size of 50 pixels excludes very small specks that are likely noise or staining artifacts rather than real signal. Adjust this value depending on the scale and resolution of your images.
>    >
>    {: .comment}
>
{: .hands-on}

This produces two outputs per image: an outline image showing detected particles, and a tabular file with ROI coordinates. Use the outline image to quickly check whether the detected regions match the staining you see in the original.

> <tip-title>Uploading ROIs to OMERO for validation</tip-title>
>
> If your facility uses OMERO for image management, the ROI files generated here can be uploaded alongside the original whole-slide images to overlay and visually validate your threshold results at scale.
>
{: .tip}

# Step 6 — Extract Quantitative Features

With a validated binary mask, you can now measure the stained area. The feature extraction tool measures properties of the labeled regions in the mask, using the grayscale stain channel for intensity information.

> <hands-on-title>Measure stained pixel area and intensity</hands-on-title>
>
> 1. {% tool [Extract image features](toolshed.g2.bx.psu.edu/repos/imgteam/2d_feature_extraction/ip_2d_feature_extraction/0.25.2+galaxy1) %} with the following parameters:
>    - {% icon param-collection %} *"Label map"*: output of **Threshold image** {% icon tool %}
>    - *"Features to compute"*: `Use the intensity image to compute additional features`
>      - {% icon param-collection %} *"Intensity image"*: output of **Extract dataset** {% icon tool %}
>      - *"Available features"*: select `label`, `area`, `area_filled`, `mean_intensity`
>
>    > <comment-title>What do these features mean?</comment-title>
>    >
>    > - **label** — the pixel class (1 = stained region)
>    > - **area** — the number of stained pixels (this is what you will use for the percentage)
>    > - **area_filled** — area including any internal holes in the stained region
>    > - **mean_intensity** — average pixel intensity within the stained region, which can complement the area measurement as a proxy for staining strength
>    >
>    {: .comment}
>
{: .hands_on}

# Step 7 — Compile Results and Calculate Percent Stained Area

You now have, for each image: the stained pixel area (from feature extraction) and the total pixel area (from image dimensions). This final section merges all per-sample results into a single table and calculates the percentage.

## Merge the Feature Results

> <hands-on-title>Collapse per-image results into one table</hands-on-title>
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} to get the first dataset from the feature collection (used to extract the header row):
>    - {% icon param-collection %} *"Input List"*: output of **Extract image features** {% icon tool %}
>    - *"How should a dataset be selected?"*: `The first dataset`
>
> 2. {% tool [Select first](Show beginning1) %} to isolate the header row:
>    - *"Select first"*: `1`
>    - {% icon param-collection %} *"from"*: output of **Extract dataset** {% icon tool %}
>
> 3. {% tool [Select last](Show tail1) %} to extract the data row from each image's feature file:
>    - *"Select last"*: `1`
>    - {% icon param-collection %} *"from"*: output of **Extract image features** {% icon tool %}
>
> 4. {% tool [Extract element identifiers](toolshed.g2.bx.psu.edu/repos/iuc/collection_element_identifiers/collection_element_identifiers/0.0.2) %} to get sample file names:
>    - {% icon param-collection %} *"Dataset collection"*: output of **Extract image features** {% icon tool %}
>
> 5. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/5.1.0) %} to merge all data rows into one file:
>    - {% icon param-collection %} *"Collection of files to collapse"*: output of **Select last** {% icon tool %}
>
> 6. {% tool [Create text file](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_text_file_with_recurring_lines/9.5+galaxy3) %} to create a `sample_id` column header:
>    - *"Line"*: `sample_id`
>
> 7. {% tool [Paste](Paste1) %} to combine the sample_id header with the feature header:
>    - {% icon param-file %} *"Paste"*: output of **Create text file** (sample_id) {% icon tool %}
>    - {% icon param-file %} *"and"*: output of **Select first** (feature header) {% icon tool %}
>
> 8. {% tool [Paste](Paste1) %} to combine sample names with their corresponding data rows:
>    - {% icon param-file %} *"Paste"*: output of **Extract element identifiers** {% icon tool %}
>    - {% icon param-file %} *"and"*: output of **Collapse Collection** {% icon tool %}
>
> 9. {% tool [Concatenate multiple datasets or collections](cat1) %} to build the full table with header:
>    - {% icon param-file %} *"Concatenate Dataset"*: output of **Paste** (header row) {% icon tool %}
>    - In *"Dataset"* → *"Insert Dataset"*:
>      - {% icon param-file %} *"Select"*: output of **Paste** (data rows) {% icon tool %}
>
{: .hands_on}

## Add the Total Area Column

> <hands-on-title>Merge the total area values</hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} to extract the total area column (column 3) from the computed area file:
>    - *"Cut columns"*: `c3`
>    - {% icon param-collection %} *"From"*: output of **Compute** (width × height) {% icon tool %}
>
> 2. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/5.1.0) %} to merge the total area values across samples:
>    - {% icon param-collection %} *"Collection of files to collapse"*: output of **Cut** {% icon tool %}
>
> 3. {% tool [Create text file](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_text_file_with_recurring_lines/9.5+galaxy3) %} to create a `total_area` column header.
>    - *"Line"*: `total_area`
>
> 4. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/9.5+galaxy3) %} to prepend the header to the total area values:
>    - {% icon param-file %} *"Datasets to concatenate"*: output of **Create text file** (total_area) {% icon tool %}
>    - In *"Dataset"* → *"Insert Dataset"*:
>      - {% icon param-file %} *"Select"*: output of **Collapse Collection** (total area) {% icon tool %}
>
> 5. {% tool [Paste](Paste1) %} to join the feature results table with the total area column:
>    - {% icon param-file %} *"Paste"*: output of **Concatenate multiple datasets or collections** {% icon tool %}
>    - {% icon param-file %} *"and"*: output of **Concatenate datasets** {% icon tool %}
>
{: .hands_on}

## Calculate Percent Stained Area

You now have a single table with all the values you need. The final step divides stained area by total area and multiplies by 100.

> <hands-on-title>Compute percent stained area</hands-on-title>
>
> 1. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: output of **Paste** (full table) {% icon tool %}
>    - *"Input has a header line with column names?"*: `Yes`
>    - In *"Expressions"* → *"Insert Expressions"*:
>      - *"Add expression"*: `c4 / c6 * 100`
>      - *"The new column name"*: `percent_area`
>    - *"If an expression cannot be computed for a row"*: `Fail the entire tool run`
>
>    > <comment-title>Understanding the formula</comment-title>
>    >
>    > `c4` = stained pixel area (from feature extraction); `c6` = total pixel area (width × height). Dividing and multiplying by 100 expresses the result as a percentage. If your column order differs, verify the column numbers by inspecting the merged table before running this step.
>    >
>    {: .comment}
>
{: .hands_on}

Your final output is a TSV table with one row per sample, containing: `sample_id`, `label`, `area`, `area_filled`, `mean_intensity`, `total_area`, and `percent_area`. This file is ready for downstream statistical analysis.

> <question-title></question-title>
>
> 1. A tissue image is 2048 × 2048 pixels. The thresholded mask contains 419,430 positive pixels. What is the percent stained area?
> 2. How would you interpret this result biologically?
>
> > <solution-title></solution-title>
> >
> > 1. Total pixels = 2048 × 2048 = 4,194,304. Percent area = (419,430 / 4,194,304) × 100 ≈ **10.0%**
> >
> > 2. <div class="IHC" markdown="1">10% CD11b-positive area indicates substantial myeloid cell infiltration in the infarcted myocardium. In the context of {% cite Rettkowski2025 %}, a reduction in this value in 4-oxo RA-treated animals would suggest the treatment suppresses immune cell mobilization and limits post-MI inflammation. Always compare against appropriate controls to account for background staining.</div><div class="MT" markdown="1">10% collagen-positive area reflects moderate fibrotic remodeling. In the context of {% cite Rettkowski2025 %}, lower values in treated animals would suggest 4-oxo RA helps preserve cardiac tissue architecture by reducing fibrosis. These results are most informative when combined with CD11b quantification, since inflammation and fibrosis are closely linked post-MI.</div>
> >
> {: .solution}
>
{: .question}

# Conclusion

You have now built and run a complete image analysis workflow that takes raw histological images and produces a quantitative, reproducible measure of staining coverage. Here is a summary of what each step contributes:

| Step | What it does | Why it matters |
|---|---|---|
| Color deconvolution | Separates mixed stain signals into individual channels | Isolates your stain from the counterstain |
| Channel extraction | Selects the channel corresponding to your stain of interest | DAB (index 1) for IHC; aniline blue (index 2) for MT |
| Image dimensions | Captures total pixel area from the original image | Provides the denominator for the percentage |
| Otsu thresholding | Converts the stain channel into a binary stained/unstained mask | Automated, consistent across images |
| Analyze particles | Generates ROIs around detected stained regions | Enables visual validation of the threshold |
| Feature extraction | Measures area and intensity of stained pixels | Produces the numerator for the percentage |
| Final table + compute | Merges results and calculates percent stained area | Delivers a ready-to-analyze summary |

The same workflow applies to both IHC and Masson's Trichrome — the only difference is the channel index you select in Step 2. By adjusting this index, the workflow can be adapted to other staining combinations as long as they are compatible with the HED color space.

<div class="IHC" markdown="1">

In the context of this tutorial, the CD11b-positive area quantified from IHC images reflects myeloid cell infiltration in infarcted cardiac tissue — an objective, reproducible measure that can be compared across experimental groups and applied consistently at scale, as demonstrated in {% cite Rettkowski2025 %}.

</div>

<div class="MT" markdown="1">

In the context of this tutorial, the collagen-positive area from Masson's Trichrome images reflects fibrotic remodeling post-MI. Combined with CD11b quantification, it provides a complementary picture of the inflammatory and structural response to infarction, as studied in {% cite Rettkowski2025 %}.

</div>

![Overview of the complete stain quantification workflow from raw image to percent stained area table.](../../images/workflow_overview.png "Workflow overview: from raw histological image to percent stained area.")
