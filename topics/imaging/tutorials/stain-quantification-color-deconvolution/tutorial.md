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
- The DAB channel is used for IHC quantification, while the aniline blue channel is used for Masson's Trichrome
- Otsu thresholding provides an automated, reproducible way to distinguish stained from unstained pixels
- The percent stained area is calculated as (stained pixels / total tissue pixels) × 100
- This workflow input is a collection of images for batch quantification
contributors:
- dianichj
---

Histological staining is widely used in biomedical research and pathology to visualize specific tissue components, proteins, or structures. However, visual assessment of staining intensity and coverage is inherently subjective and difficult to reproduce across studies. Computational quantification of stained area provides an objective, reproducible measurement that can be applied consistently across large image datasets and used to address different research questions.

**Color deconvolution** is a key technique for this purpose. Originally described by {% cite Ruifrok2001 %}, it mathematically separates the individual stain components in a brightfield RGB image based on their characteristic optical absorption spectra (stain vectors). This allows you to isolate one stain from the counterstain and independently quantify its coverage. For example, in Immunohistochemistry (IHC) it is possible to separate the DAB (3,3'-diaminobenzidine) signal from the hematoxylin counterstain.

In this tutorial, you will learn how to use Galaxy's image analysis tools to:

1. Apply color deconvolution to histological images using predefined stain presets
2. Isolate the channel of interest for IHC (DAB) or Masson's Trichrome (aniline blue)
3. Automatically threshold the deconvolved channel to detect positive staining
4. Extract image features and calculate the percentage of stained area
5. Compile results across multiple images into a final summary table

This workflow is applicable to **IHC staining** (e.g. antibody detection with DAB chromogen) and **Masson's Trichrome (MT)** staining (e.g. collagen quantification), and can be extended to other staining combinations by changing the deconvolution preset and adjusting the threshold type accordingly. It was originally inspired by the ImageJ-based collagen quantification method described by {% cite Chen2017 %}, adapted for cardiac tissue section analysis, and subsequently implemented in Galaxy to enable reproducible, scalable batch processing.

> <agenda-title></agenda-title>
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Background: Color Deconvolution

Before starting the hands-on part, it is important to understand what color deconvolution does and why it is needed.

In brightfield histology, stains absorb light at specific wavelengths, producing characteristic colors in the RGB image. When two or more stains are applied simultaneously (e.g. hematoxylin + DAB, or Weigert's iron hematoxylin + Biebrich scarlet + aniline blue), their colors overlap in the RGB channels. Simple color thresholding on the raw RGB image cannot reliably separate them.

Color deconvolution solves this by using a **stain matrix**, a set of stain vectors that describe the optical density contributions of each stain across the three RGB channels. The algorithm inverts this matrix to compute the individual optical density of each stain at every pixel, producing a separate grayscale image per stain channel.

**Common stain presets available in Galaxy:**

| Preset | Stain 1 | Stain 2 | Stain 3 | Use case |
|---|---|---|---|---|
| H + DAB | Hematoxylin | DAB | (residual) | IHC |
| H + E + DAB | Hematoxylin | Eosin | DAB | IHC with H&E counterstain |
| Masson's Trichrome | Weigert's hematoxylin | Biebrich scarlet | Aniline blue | Collagen quantification |

> <details-title> What is a stain vector? </details-title>
>
> A stain vector describes how much a given stain absorbs light in each of the three RGB color channels. For example, DAB absorbs strongly in the blue channel and weakly in the red channel, while hematoxylin absorbs more in the red and green channels. These values are determined empirically from pure stain reference images or taken from published standards. The color deconvolution algorithm uses these vectors to solve a system of linear equations and separate the mixed stain signals pixel by pixel.
>
{: .details}

# Get Data

The images used in this tutorial are a subset of data from {% cite Rettkowski2025 %}, a study investigating the role of 4-oxo retinoic acid (4-oxo RA) in maintaining hematopoietic stem cell dormancy in the bone marrow following myocardial infarction (MI). The biological model consists of mouse hearts subjected to left anterior descending (LAD) coronary artery ligation to induce MI, with and without 4-oxo RA treatment.

Two staining types were applied to serial cardiac sections:

- **IHC for CD11b** — a surface marker of myeloid leukocytes (monocytes, macrophages, neutrophils) used here as a proxy for local inflammation. In this context, higher CD11b-positive area is associated with greater immune cell presence post-MI.
- **Masson's Trichrome (MT)** — to visualize and quantify collagen deposition, a hallmark of cardiac fibrosis. Greater collagen content reflects more extensive fibrotic remodeling following infarction.

The hypothesis is that 4-oxo RA treatment reduces immune cell mobilization from the bone marrow, leading to less CD11b-positive infiltration, less fibrosis, and better preservation of cardiac tissue post-MI. Images were acquired from three anatomical zones of the heart (Infarct, Border, and Remote), although zone is not used as a variable in this tutorial.

In the original study, regions of interest (ROIs) were manually selected from whole-slide images to focus the analysis on specific anatomical zones and avoid tissue artifacts. For this tutorial, the ROIs have already been cropped and are provided directly as the input image collection, so you can focus on the quantification workflow itself.

> <comment-title> Sample data vs. full dataset </comment-title>
>
> The images provided here are a representative sample of the full dataset used in {% cite Rettkowski2025 %}. In practice, histological image quantification workflows like this one are run on large batches of images across multiple animals, conditions, and zones. The workflow you will run here is identical to what would be applied at scale.
>
{: .comment}

## Data upload

> <hands-on-title> Data Upload </hands-on-title>
>
> 1. Create a new history for this tutorial and name it something meaningful, e.g. `Histology Stain Quantification`
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
> 4. Check that the datatype is set to `tiff` or `png`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tiff" %}
>
> 5. Organize your images into a **dataset collection** (one collection per staining type). This allows the workflow to process all images in batch.
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
> > <comment-title> Image requirements </comment-title>
> >
> > Images must be **brightfield RGB** microscopy images. Fluorescence images are not compatible with color deconvolution. Avoid images with strong artifacts, out-of-focus regions, or uneven illumination, as these will reduce quantification accuracy.
> >
> {: .comment}
>
{: .hands_on}

{% include _includes/cyoa-choices.html option1="IHC" option2="MT" default="IHC"
   text="This tutorial can be followed with two different staining types: Immunohistochemistry (IHC) for CD11b, or Masson's Trichrome (MT) for collagen quantification. Choose the one that matches your data!" %}

# Color Deconvolution and Channel Extraction

In this section you will apply color deconvolution to your input images and isolate the channel corresponding to your stain of interest.

<div class="IHC" markdown="1">

In IHC images stained with DAB chromogen, color deconvolution separates the brown DAB signal from the blue hematoxylin counterstain, producing one grayscale channel per stain component. You will use the DAB channel to quantify CD11b-positive area.

> <hands-on-title> Perform Color Deconvolution (IHC) </hands-on-title>
>
> 1. {% tool [Perform color deconvolution or transformation](toolshed.g2.bx.psu.edu/repos/imgteam/color_deconvolution/ip_color_deconvolution/0.9+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input image"*: your IHC image collection
>    - *"Transformation type"*: `Deconvolve RGB into Hematoxylin + DAB`
>
>    > <comment-title> Output </comment-title>
>    >
>    > The tool produces one multi-channel TIFF image per input image. Each channel is a grayscale image representing the optical density of one stain component: brighter pixels indicate stronger staining. Channel 1 corresponds to Hematoxylin and channel 2 to DAB.
>    >
>    {: .comment}
>
{: .hands_on}

![Color deconvolution output for IHC (CD11b/DAB). From left to right: original RGB image, Hematoxylin channel (channel 1), and DAB channel (channel 2).](../../images/deconvolution_output_IHC.png "Color deconvolution output for IHC. Left: original RGB image. Middle: Hematoxylin channel. Right: DAB channel.")

> <details-title> When to use Non-negative Matrix Factorization (NMF) instead </details-title>
>
> For most standard IHC images with a clean DAB + hematoxylin combination, the `Deconvolve RGB into Hematoxylin + DAB` preset is sufficient. However, in cases where staining is uneven, the DAB signal is weak, or there is significant spectral overlap between stain components, a data-driven approach may give better results.
>
> **Non-negative Matrix Factorization (NMF)** is an unsupervised method that learns the stain components directly from the image data, without relying on predefined stain vectors. Instead of using fixed optical absorption spectra, NMF decomposes the image into a set of non-negative components that best explain the observed pixel values. This makes it more flexible and potentially more accurate when the staining deviates from the standard reference spectra — for example, due to differences in staining protocols, tissue processing, or scanner calibration across batches.
>
> To use NMF, select `Non-negative matrix factorization` as the transformation type. Note that since NMF components are not labeled by stain name, you will need to visually inspect the output channels to identify which one corresponds to your stain of interest before proceeding.
>
{: .details}

</div>

<div class="MT" markdown="1">

In Masson's Trichrome images, color deconvolution separates the blue aniline blue signal (collagen) from the red Biebrich scarlet (muscle/cytoplasm) and the dark Weigert's hematoxylin (nuclei). You will use the third channel to quantify collagen-positive area.

> <hands-on-title> Perform Color Deconvolution (MT) </hands-on-title>
>
> 1. {% tool [Perform color deconvolution or transformation](toolshed.g2.bx.psu.edu/repos/imgteam/color_deconvolution/ip_color_deconvolution/0.9+galaxy0) %} with the following parameters:
>    - {% icon param-collection %} *"Input image"*: your MT image collection
>    - *"Transformation type"*: `Deconvolve RGB into Hematoxylin + Eosin + DAB`
>
>    > <comment-title> Output </comment-title>
>    >
>    > The tool produces one multi-channel TIFF image per input image. Although this preset is named after H&E+DAB, its stain vectors provide a good separation for Masson's Trichrome components. Channel 1 corresponds to Hematoxylin, channel 2 to Eosin/Biebrich scarlet, and channel 3 to the collagen/aniline blue signal.
>    >
>    {: .comment}
>
{: .hands_on}

![Color deconvolution output for Masson's Trichrome. From left to right: original RGB image, Hematoxylin channel (channel 1), Eosin/Biebrich scarlet channel (channel 2), and collagen/aniline blue channel (channel 3).](../../images/deconvolution_output_MT.png "Color deconvolution output for Masson's Trichrome. Left: original RGB image. Channels 1-3 shown separately.")

</div>

# Image Splitting and Channel Extraction

Once color deconvolution is complete, you will split the multi-channel output into individual images and extract the channel of interest.

> <hands-on-title> Split image along axes </hands-on-title>
>
> 1. {% tool [Split image along axes](toolshed.g2.bx.psu.edu/repos/imgteam/split_image/ip_split_image/2.3.5+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Image to split"*: `output` (output of **Perform color deconvolution or transformation** {% icon tool %})
>
>    > <comment-title> Output </comment-title>
>    >
>    > This tool splits the multi-channel TIFF into a collection of single-channel grayscale images, one per stain component. You will select the channel of interest in the next step.
>    >
>    {: .comment}
>
{: .hands_on}

<div class="IHC" markdown="1">

> <hands-on-title> Extract DAB channel (IHC) </hands-on-title>
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-file %} *"Input List"*: `output` (output of **Split image along axes** {% icon tool %})
>    - *"How should a dataset be selected?"*: `Select by index`
>    - *"Element index"*: `1`
>
>    > <comment-title> Why index 1? </comment-title>
>    >
>    > The split collection is zero-indexed. In the H + DAB preset, DAB is the second channel, corresponding to index `1`.
>    >
>    {: .comment}
>
{: .hands_on}

</div>

<div class="MT" markdown="1">

> <hands-on-title> Extract collagen channel (MT) </hands-on-title>
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} with the following parameters:
>    - {% icon param-file %} *"Input List"*: `output` (output of **Split image along axes** {% icon tool %})
>    - *"How should a dataset be selected?"*: `Select by index`
>    - *"Element index"*: `2`
>
>    > <comment-title> Why index 2? </comment-title>
>    >
>    > The split collection is zero-indexed. In the H + E + DAB preset used for MT, the collagen/aniline blue signal is the third channel, corresponding to index `2`.
>    >
>    {: .comment}
>
{: .hands_on}

</div>

> <question-title></question-title>
>
> 1. What does a high pixel intensity value mean in the deconvolved DAB channel?
> 2. What index would you use to extract the hematoxylin channel from an H + DAB deconvolution?
>
> > <solution-title></solution-title>
> >
> > 1. After deconvolution, the output is an optical density image. Higher pixel values represent stronger DAB staining, meaning more DAB chromogen is present at that location. Darker regions in the original image correspond to higher values in the deconvolved channel.
> > 2. Hematoxylin is the first channel, so you would use index `0`.
> >
> {: .solution}
>
{: .question}

# Retrieve Image Dimensions

To calculate the percentage of stained area, you need two numbers: the area covered by the stain (which you will measure later from the thresholded image), and the total area of the image, which serves as the denominator. In this section you will extract the width and height of each image from its metadata and multiply them to obtain the total pixel area.

> <hands-on-title> Get image dimensions </hands-on-title>
>
> 1. {% tool [Show image info](toolshed.g2.bx.psu.edu/repos/imgteam/image_info/ip_imageinfo/5.7.1+galaxy1) %} with the following parameters:
>    - {% icon param-collection %} *"Input Image"*: `output` (output of **Perform color deconvolution or transformation** {% icon tool %})
>
> 2. {% tool [Select](Grep1) %} to extract the image width:
>    - {% icon param-file %} *"Select lines from"*: `output` (output of **Show image info** {% icon tool %})
>    - *"the pattern"*: `Width =`
>
> 3. {% tool [Select](Grep1) %} to extract the image height:
>    - {% icon param-file %} *"Select lines from"*: `output` (output of **Show image info** {% icon tool %})
>    - *"the pattern"*: `Height =`
>
> 4. {% tool [Text transformation](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sed_tool/9.5+galaxy3) %} to isolate the width value:
>    - {% icon param-file %} *"File to process"*: output of **Select** (Width) {% icon tool %}
>    - *"SED Program"*: `s/.*= //`
>    - *"Advanced Options"*: `Hide Advanced Options`
>
>    > <comment-title> What does this SED expression do? </comment-title>
>    >
>    > The **Show image info** tool returns metadata as plain text lines, e.g. `Width = 1024`. The expression `s/.*= //` removes everything up to and including the `= ` sign, leaving only the numeric value. This is necessary so the number can be used in a calculation in the next steps.
>    >
>    {: .comment}
>
> 5. {% tool [Text transformation](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_sed_tool/9.5+galaxy3) %} to isolate the height value:
>    - {% icon param-file %} *"File to process"*: output of **Select** (Height) {% icon tool %}
>    - *"SED Program"*: `s/.*= //`
>    - *"Advanced Options"*: `Hide Advanced Options`
>
> 6. {% tool [Paste](Paste1) %} to combine width and height into a single file:
>    - {% icon param-file %} *"Paste"*: output of **Text transformation** (width) {% icon tool %}
>    - {% icon param-file %} *"and"*: output of **Text transformation** (height) {% icon tool %}
>
> 7. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.1+galaxy0) %} to calculate total pixel area (width × height):
>    - {% icon param-file %} *"Input file"*: output of **Paste** {% icon tool %}
>    - *"Input has a header line with column names?"*: `No`
>    - In *"Expressions"* → *"Insert Expressions"*:
>      - *"Add expression"*: `c1 * c2`
>    - *"If an expression cannot be computed for a row"*: `Fail the entire tool run`
>
{: .hands_on}

> <question-title></question-title>
>
> An image is 1024 × 768 pixels. What is the total pixel area, and why does it matter for this analysis?
>
> > <solution-title></solution-title>
> >
> > The total pixel area is 1024 × 768 = 786,432 pixels. This value is the denominator in the final percentage calculation: stained pixels / total pixels × 100. Without it, you can measure the absolute number of stained pixels but cannot express it as a meaningful percentage that is comparable across images of different sizes.
> >
> {: .solution}
>
{: .question}

# Thresholding and Feature Extraction

Once you have isolated the stain channel, you need to convert it into a binary mask (stained vs. unstained pixels) and then measure the area covered by the stain.

## Threshold the Stain Channel

> <hands-on-title> Apply automatic thresholding </hands-on-title>
>
> 1. {% tool [Threshold image](toolshed.g2.bx.psu.edu/repos/imgteam/2d_auto_threshold/ip_threshold/0.25.2+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input image"*: `output` (output of **Extract dataset** {% icon tool %})
>    - *"Thresholding method"*: `Globally adaptive / Otsu`
>
>    > <comment-title> Why Otsu thresholding? </comment-title>
>    >
>    > Otsu's method automatically finds an optimal threshold by minimizing the intra-class variance between the two pixel populations (stained and unstained). This makes it more reproducible than manual threshold selection, as it does not require user input per image. It works best when the pixel intensity histogram shows two distinct peaks (bimodal distribution), which is typically the case for well-stained histological sections.
>    >
>    {: .comment}
>
{: .hands_on}

The figure below shows an example of the thresholding output. The binary mask highlights the stained pixels in white against a black background:

<div class="IHC" markdown="1">

![Thresholding output for IHC. Left: deconvolved DAB channel. Right: binary mask after Otsu thresholding, where white pixels represent DAB-positive area.](../../images/threshold_output_IHC.png "Thresholding output for IHC (CD11b/DAB). Left: DAB channel. Right: Otsu binary mask.")

</div>

<div class="MT" markdown="1">

![Thresholding output for Masson's Trichrome. Left: deconvolved collagen channel. Right: binary mask after Otsu thresholding, where white pixels represent collagen-positive area.](../../images/threshold_output_MT.png "Thresholding output for MT. Left: collagen channel. Right: Otsu binary mask.")

</div>

> <question-title></question-title>
>
> 1. What does the output of the threshold step look like?
> 2. When might Otsu thresholding give unreliable results?
>
> > <solution-title></solution-title>
> >
> > 1. The output is a binary mask image: pixels classified as stained are set to `1` (white) and unstained pixels are set to `0` (black).
> > 2. Otsu thresholding may underperform when the image has very low or very high staining density (unimodal histogram), or when there is significant background noise or uneven illumination. In those cases, consider preprocessing images (e.g. background subtraction) or using a manually set threshold.
> >
> {: .solution}
>
{: .question}

## Extract Image Features from the Binary Mask

> <hands-on-title> Measure stained pixel area </hands-on-title>
>
> 1. {% tool [Extract image features](toolshed.g2.bx.psu.edu/repos/imgteam/2d_feature_extraction/ip_2d_feature_extraction/0.25.2+galaxy1) %} with the following parameters:
>    - {% icon param-file %} *"Label map"*: `output` (output of **Threshold image** {% icon tool %})
>    - *"Features to compute"*: `Use the intensity image to compute additional features`
>      - {% icon param-file %} *"Intensity image"*: `output` (output of **Extract dataset** {% icon tool %})
>      - *"Available features"*: *(select area and intensity features)*
>
>    > <comment-title> What features are extracted? </comment-title>
>    >
>    > This tool computes morphological and intensity features from labeled regions in the binary mask. The key output here is **area** — the number of pixels classified as stained. This value will be used in the final calculation of percentage stained area.
>    >
>    {: .comment}
>
{: .hands_on}

# Compiling and Calculating Percentage Stained Area

In this final section, you will gather the stained area values and total image area across all images and compute the percentage stained area for each sample.

## Collect and Organize Feature Data

> <hands-on-title> Extract and collapse feature values </hands-on-title>
>
> 1. {% tool [Extract dataset](__EXTRACT_DATASET__) %} to get the first dataset from the feature collection:
>    - {% icon param-file %} *"Input List"*: `output` (output of **Extract image features** {% icon tool %})
>    - *"How should a dataset be selected?"*: `The first dataset`
>
> 2. {% tool [Select last](Show tail1) %} to get the last row (the stained region values):
>    - *"Select last"*: `1`
>    - {% icon param-file %} *"from"*: `output` (output of **Extract image features** {% icon tool %})
>
> 3. {% tool [Extract element identifiers](toolshed.g2.bx.psu.edu/repos/iuc/collection_element_identifiers/collection_element_identifiers/0.0.2) %} to retrieve image file names:
>    - {% icon param-file %} *"Dataset collection"*: `output` (output of **Extract image features** {% icon tool %})
>
> 4. {% tool [Cut](Cut1) %} to extract the stained area column:
>    - *"Cut columns"*: `c3`
>    - {% icon param-file %} *"From"*: `out_file1` (output of **Compute** {% icon tool %})
>
> 5. {% tool [Select first](Show beginning1) %} to get the header row:
>    - *"Select first"*: `1`
>    - {% icon param-file %} *"from"*: `output` (output of **Extract dataset** {% icon tool %})
>
{: .hands_on}

## Collapse Collections into Summary Tables

> <hands-on-title> Merge results across all images </hands-on-title>
>
> 1. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/5.1.0) %} to merge stained area rows:
>    - {% icon param-file %} *"Collection of files to collapse into single dataset"*: `out_file1` (output of **Select last** {% icon tool %})
>
> 2. {% tool [Collapse Collection](toolshed.g2.bx.psu.edu/repos/nml/collapse_collections/collapse_dataset/5.1.0) %} to merge total area column:
>    - {% icon param-file %} *"Collection of files to collapse into single dataset"*: `out_file1` (output of **Cut** {% icon tool %})
>
> 3. {% tool [Paste](Paste1) %} to join image names with stained area:
>    - {% icon param-file %} *"Paste"*: `output` (output of **Extract element identifiers** {% icon tool %})
>    - {% icon param-file %} *"and"*: `output` (output of **Collapse Collection** — stained area {% icon tool %})
>
> 4. {% tool [Paste](Paste1) %} to add a header to the table:
>    - {% icon param-file %} *"Paste"*: `outfile` (output of **Create text file** {% icon tool %})
>    - {% icon param-file %} *"and"*: `out_file1` (output of **Select first** {% icon tool %})
>
> 5. {% tool [Concatenate datasets](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/9.5+galaxy3) %} to build the total area column:
>    - {% icon param-file %} *"Datasets to concatenate"*: `outfile` (output of **Create text file** {% icon tool %})
>    - In *"Dataset"* → *"Insert Dataset"*:
>      - {% icon param-file %} *"Select"*: `output` (output of **Collapse Collection** — total area {% icon tool %})
>
> 6. {% tool [Concatenate multiple datasets or collections](cat1) %} to assemble the full table:
>    - {% icon param-file %} *"Concatenate Dataset"*: `out_file1` (output of **Paste** — names + stained area {% icon tool %})
>    - In *"Dataset"* → *"Insert Dataset"*:
>      - {% icon param-file %} *"Select"*: `out_file1` (output of **Paste** — header {% icon tool %})
>
> 7. {% tool [Paste](Paste1) %} to combine all columns:
>    - {% icon param-file %} *"Paste"*: `out_file1` (output of **Concatenate multiple datasets or collections** {% icon tool %})
>    - {% icon param-file %} *"and"*: `out_file1` (output of **Concatenate datasets** {% icon tool %})
>
{: .hands_on}

## Calculate Percentage Stained Area

This is the final step: dividing the stained pixel area by the total image area and multiplying by 100 to obtain the percentage.

> <hands-on-title> Compute percent stained area </hands-on-title>
>
> 1. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.1+galaxy0) %} with the following parameters:
>    - {% icon param-file %} *"Input file"*: `out_file1` (output of **Paste** {% icon tool %})
>    - *"Input has a header line with column names?"*: `Yes`
>    - In *"Expressions"* → *"Insert Expressions"*:
>      - *"Add expression"*: `c4 / c6 * 100`
>      - *"The new column name"*: `percent_area`
>    - *"If an expression cannot be computed for a row"*: `Fail the entire tool run`
>
>    > <comment-title> Understanding the formula </comment-title>
>    >
>    > The expression `c4 / c6 * 100` divides the stained pixel area (column 4) by the total image area (column 6) and multiplies by 100 to express the result as a percentage. Make sure your column assignments match the actual structure of your merged table.
>    >
>    {: .comment}
>
{: .hands_on}

> <question-title></question-title>
>
> 1. A tissue image is 2048 × 2048 pixels. The thresholded mask contains 419,430 positive pixels. What is the percent stained area?
> 2. How would you interpret the results biologically?
>
> > <solution-title></solution-title>
> >
> > 1. Total pixels = 2048 × 2048 = 4,194,304. Percent area = (419,430 / 4,194,304) × 100 ≈ 10.0%
> >
> > 2. <div class="IHC" markdown="1">A higher CD11b-positive area indicates greater myeloid cell infiltration, which reflects local inflammation in the infarcted myocardium. In the context of {% cite Rettkowski2025 %}, a reduction in CD11b-positive area in 4-oxo RA-treated animals would suggest that the treatment reduces immune cell mobilization from the bone marrow and limits inflammatory infiltration post-MI. Always compare against negative controls to account for non-specific background staining.</div><div class="MT" markdown="1">A higher collagen-positive area reflects greater fibrotic remodeling of the myocardium following infarction. In the context of {% cite Rettkowski2025 %}, reduced collagen deposition in treated animals would suggest better preservation of cardiac tissue architecture. Results should be interpreted alongside CD11b quantification, as reduced inflammation is expected to correlate with reduced fibrosis.</div>
> >
> {: .solution}
>
{: .question}

# Conclusion

In this tutorial, you have learned how to use Galaxy's image analysis tools to quantitatively measure the percentage of stained area in histological images. The key steps of the workflow are:

1. **Color deconvolution** — separates mixed stain signals into individual channels using optical absorption spectra
2. **Channel extraction** — isolates the stain of interest (DAB for IHC, aniline blue for MT)
3. **Otsu thresholding** — automatically creates a binary mask of stained vs. unstained pixels
4. **Feature extraction** — measures the pixel area covered by the stain
5. **Percentage calculation** — expresses the stained area relative to total image area as a percentage

This approach provides an objective, reproducible, and scalable method for histological stain quantification that can be applied to large image datasets. By adjusting the color deconvolution preset and channel index, the same workflow can be adapted to other staining protocols beyond IHC and Masson's Trichrome.

<div class="IHC" markdown="1">

In the context of this tutorial, the CD11b-positive area quantified from IHC images reflects the degree of myeloid cell infiltration in infarcted cardiac tissue. This type of measurement allows objective comparison between experimental groups — for example, treated vs. untreated animals — and can be applied consistently across large image collections from studies such as {% cite Rettkowski2025 %}.

</div>

<div class="MT" markdown="1">

In the context of this tutorial, the collagen-positive area quantified from Masson's Trichrome images reflects the degree of fibrotic remodeling in infarcted cardiac tissue. Combined with CD11b quantification, this measurement provides a complementary readout of the inflammatory and fibrotic response post-MI, as studied in {% cite Rettkowski2025 %}.

</div>

![Workflow overview: from raw histological image to percent stained area.](../../images/workflow_overview.png "Overview of the stain quantification workflow.")