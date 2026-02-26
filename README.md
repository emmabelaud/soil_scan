# In situ scanner-based imaging to monitor soil edaphic activity

This repository supports the development of a soil biological activity monitoring method based on **buried flatbed scanners**. The approach relies on CIS-type flatbed scanners installed vertically in the soil, which continuously capture high-resolution images of the soil matrix *in situ*. These images enable the observation and tracking of various soil organisms — including invertebrates and roots structures — at fine spatial and temporal resolutions. Coupled with automated image capture (via Raspberry Pi modules) and a suite of image processing and analysis tools, this method offers a non-destructive, long-term window into belowground biological dynamics that traditional sampling approaches cannot easily provide.

The methodology spans two main components: **image acquisition** (scanner setup, automated capture, and parameterization) and **image analysis** (detection and classification of invertebrates and roots, using both manual workflows and automated deep learning tools).

This non-invasive, temporally continuous approach opens new possibilities for studying edaphic biodiversity across diverse spatio-temporal scales ([Belaud et al. 2024](https://doi.org/10.1007/s00374-024-01851-8)).

![](image/scan-device.png)

## Image acquisition

The [image_acquisition](https://github.com/emmabelaud/soil_scan/tree/30b0c6e3ef320dadd86181c5e8fc15f40663755e/images_acquisition) folder contains all code and data associated with the paper *"In situ scanner-based imaging to monitor soil mesofauna activity: trade-offs between spatial resolution and temporal frequency"*.

A core challenge in deploying buried scanner systems for ecological monitoring is that image-derived activity metrics are not straightforward proxies of biological activity — they are inherently shaped by how images are acquired. This study addresses this critical but rarely quantified dimension by experimentally **calibrating the effects of two key acquisition parameters — spatial resolution and capture frequency — on the detectability of soil mesofauna** under controlled greenhouse conditions.

Focusing on three major taxonomic groups (Acari, Collembola, and Enchytraeidae), it demonstrates that spatial resolution acts as the dominant driver of detectability, functioning as a size-selective filter that determines which organisms and size classes enter the observation window. Acquisition frequency, in turn, modulates detection magnitude through temporal integration without redefining the detectable size spectrum. Together, these findings provide a mechanistic framework for understanding how protocol choices translate into observable signals, and offer practical guidance for standardising acquisition settings across studies — a prerequisite for robust ecological inference and cross-study comparability in automated soil biodiversity monitoring.

## Image analysis

The [image_analysis](https://github.com/emmabelaud/soil_scan/tree/30b0c6e3ef320dadd86181c5e8fc15f40663755e/images_analysis) folder contains all code and data associated with the paper *"Automated processing for in situ soil monitoring: from raw images to population dynamics"*. [View the analysis.](https://htmlpreview.github.io/?https://github.com/emmabelaud/soil_scan/blob/24b6554309cb9a2cf67c79a98833a7e2e406e59e/images_analysis/Use-case.html)

While buried scanner systems can generate image datasets of unprecedented spatio-temporal depth, translating this raw visual data into exploitable biological information represents a significant processing bottleneck. This study tackles this challenge by developing an **efficient end-to-end processing pipeline that combines traditional computer vision algorithms with state-of-the-art deep learning models**, specifically designed to operate under the constrained annotation budgets typical of ecological studies.

A key methodological contribution is the use of image differencing to reframe the initial detection problem — complicated by very low signal-to-noise ratios in soil imagery — into a simpler classification task, solved by fine-tuning foundation models on limited labelled data. This approach enabled the release of a dataset comprising approximately 600 soil images and over 8 000 labelled invertebrates, spanning taxa with highly varying shapes and sizes (Acari, Collembola, Enchytraeidae, and others). The resulting pipeline delivers reliable population counts across taxa and was applied to analyse soil fauna population dynamics across four scanners over 3 month of continuous monitoring — demonstrating that thoughtfully combining heuristic domain knowledge with modern AI methods can substantially reduce annotation burden while unlocking longitudinal biodiversity analyses at scales previously out of reach.
