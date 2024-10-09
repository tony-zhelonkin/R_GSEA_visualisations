# Custom Pathway Visualizations Function

This repository contains custom R functions for comprehensive pathway analysis visualisation. They are inteded to be used on the GSEA, gseGO and other pathway analysis methods results from the `ClusterProfiler` R library.

The functions are:
* custom_dotplot
* custom_dotplot_comparison
* gsea_comparison_plot
* analyze_pathway_volcano_custom

## Custom GSEA Dotplot Function

This repository contains an R function `custom_dotplot` for creating customized dotplots from Gene Set Enrichment Analysis (GSEA) results.

### Function Description

The `custom_dotplot` function generates a visually appealing and informative dotplot from GSEA results, allowing for various customizations including text processing, filtering, and sorting options.

#### Parameters

- `gsea_obj`: A GSEA result object
- `showCategory`: Number of categories to display (default: 10)
- `font.size`: Base font size for the plot (default: 10)
- `title`: Plot title (default: "GSEA Dotplot")
- `replace_`: Boolean to replace underscores with spaces in descriptions (default: TRUE)
- `capitalize_1`: Boolean to capitalize the first word in descriptions (default: TRUE)
- `capitalize_all`: Boolean to capitalize all words in descriptions (default: FALSE)
- `filterBy`: Method to filter results ("qvalue", "NES", "NES_positive", "NES_negative")
- `sortBy`: Method to sort results ("GeneRatio" or "qvalue")
- `q_cut`: Q-value cutoff for significance (default: 0.05)
- `min.dotSize`: Minimum dot size in the plot (default: 2)

#### Dependencies

The function requires the following R packages:

- ggplot2
- dplyr
- stringr
- scales

### Usage

```R
source("custom_dotplot.R")

# Assuming you have a GSEA result object named 'gsea_result'
plot <- custom_dotplot(gsea_result, showCategory = 15, filterBy = "NES", sortBy = "GeneRatio")
print(plot)
:```

### Features

1. Flexible filtering options based on q-value, NES (Normalized Enrichment Score), or direction of enrichment
2. Customizable sorting of results by Gene Ratio or q-value
3. Text processing options for pathway descriptions (capitalization, underscore replacement)
4. Visually appealing dotplot with:
   - Dot size representing significance (-log10(q-value))
   - Color indicating direction of enrichment (positive or negative NES)
   - Customizable appearance (font size, title, etc.)

### Customization

You can customize the plot by adjusting various parameters:

- `showCategory` to control the number of pathways displayed
- `filterBy` and `sortBy` for different result filtering and sorting strategies
- `q_cut` to set the significance threshold
- `font.size` and `title` for visual customization

### Contributing

Feel free to fork this repository and submit pull requests with improvements or bug fixes.

### License

[Specify your chosen license here]




## Custom GSEA Dotplot Comparison Function

This repository contains an R function `custom_dotplot_comparison` for creating a comparative dotplot visualization of Gene Set Enrichment Analysis (GSEA) results from two samples.

### Function Description

The `custom_dotplot_comparison` function generates a side-by-side dotplot comparison of GSEA results from two samples, allowing for easy visualization of enriched pathways across both datasets.

#### Parameters

- `gsea_obj_x`, `gsea_obj_y`: GSEA result objects for comparison
- `pathway_ids`: Vector of pathway IDs to include in the plot
- `font.size`: Base font size for the plot (default: 10)
- `title`: Plot title (default: "Combined GSEA Dotplot")
- `replace_`: Boolean to replace underscores with spaces in descriptions (default: TRUE)
- `capitalize_1`: Boolean to capitalize the first word in descriptions (default: FALSE)
- `capitalize_all`: Boolean to capitalize all words in descriptions (default: FALSE)
- `min.dotSize`: Minimum dot size in the plot (default: 2)
- `sample_x_name`, `sample_y_name`: Names for the two samples (default: "Sample X" and "Sample Y")
- `sortBy`: Method to sort pathways ("importance_score", "qvalue", "x", or "y") (default: "importance_score")

#### Dependencies

The function requires the following R packages:

- ggplot2
- dplyr
- stringr
- tidyr

### Usage

```R
source("custom_dotplot_comparison.R")

# Assuming you have two GSEA result objects named 'gsea_result1' and 'gsea_result2'
# and a vector of pathway IDs you want to compare
plot <- custom_dotplot_comparison(
  gsea_obj_x = gsea_result1,
  gsea_obj_y = gsea_result2,
  pathway_ids = c("PATHWAY1", "PATHWAY2", "PATHWAY3"),
  sample_x_name = "Treatment",
  sample_y_name = "Control",
  sortBy = "importance_score"
)
print(plot)
```

### Features

1. Generates a side-by-side dotplot comparison of GSEA results from two samples
2. Allows customization of pathway descriptions (replacing underscores, capitalization)
3. Provides flexible sorting options for pathways
4. Uses color to distinguish between positive and negative Normalized Enrichment Scores (NES)
5. Represents significance with dot size (-log10(q-value))
6. Handles missing pathways gracefully with warnings

### Customization

You can customize the plot by adjusting various parameters:

- `font.size` for changing the base font size
- `replace_`, `capitalize_1`, and `capitalize_all` for pathway description formatting
- `min.dotSize` to set the minimum dot size
- `sortBy` to change the ordering of pathways in the plot

### Output

The function returns a ggplot object representing the comparative dotplot. The plot includes:

- Two facets, one for each sample
- Dots representing pathways, with size indicating significance and color indicating NES direction
- Customizable pathway descriptions on the y-axis
- Gene Ratio on the x-axis

### Contributing

Feel free to fork this repository and submit pull requests with improvements or bug fixes.

### License

[Specify your chosen license here]




## GSEA NESscatter Comparison Plot Function

This repository contains an R function `gsea_comparison_plot` for comparing and visualizing Gene Set Enrichment Analysis (GSEA) results between two datasets.

### Function Description

The `gsea_comparison_plot` function generates a comprehensive comparison of GSEA results from two datasets, including a scatter plot visualization, data preparation, and extraction of common pathways.

#### Parameters

- `x`, `y`: GSEA result dataframes for comparison
- `x_label`, `y_label`: Labels for the two datasets
- `use_normalized`: Boolean to use Normalized Enrichment Score (NES) instead of Enrichment Score (default: TRUE)
- `percentile_threshold`: Threshold for determining significant pathways (default: 0.95)
- `color`: Color scheme for the plot ("all", "common", or "distinct")
- `max_overlaps`: Maximum number of label overlaps in the plot (default: 20)

#### Dependencies

The function requires the following R packages:

- dplyr
- ggplot2
- ggrepel
- tidyr

### Usage

```R
source("gsea_comparison_plot.R")

# Assuming you have two GSEA result dataframes named 'gsea_result1' and 'gsea_result2'
comparison_results <- gsea_comparison_plot(gsea_result1, gsea_result2,
                                           x_label = "Dataset1", y_label = "Dataset2",
                                           use_normalized = TRUE,
                                           percentile_threshold = 0.95,
                                           color = "all")

# Access the plot
print(comparison_results$plot)

# Access the results dataframe
head(comparison_results$data)

# Access common pathways
common_pathways <- comparison_results$common_pathways

# Access GMT-formatted dataframes
gmt_data <- comparison_results$gmt
```

### Features

1. Generates a scatter plot comparing enrichment scores or NES between two datasets
2. Identifies and categorizes common and distinct pathways
3. Provides flexible color schemes for visualizing results
4. Implements a colorblind-friendly palette
5. Labels top pathways using ggrepel for improved readability
6. Generates GMT-formatted dataframes for further analysis
7. Returns a comprehensive list object with plot, data, common pathways, and GMT-formatted results

### Customization

You can customize the plot and analysis by adjusting various parameters:

- `use_normalized` to toggle between Enrichment Score and Normalized Enrichment Score
- `percentile_threshold` to adjust the significance threshold
- `color` to change the color scheme of the plot
- `max_overlaps` to control the density of labels in the plot

### Output

The function returns a list containing:

1. `data`: A dataframe with combined results from both datasets
2. `plot`: The ggplot object for the comparison plot
3. `common_pathways`: A list of common pathways (positive, negative, mixed, and all)
4. `gmt`: A list of GMT-formatted dataframes for different pathway categories

### Contributing

Feel free to fork this repository and submit pull requests with improvements or bug fixes.

### License

[Specify your chosen license here]




## Custom Pathway Volcano Plot Function

This repository contains an R function `analyze_pathway_volcano_custom` for creating customized volcano plots for pathway analysis in the context of differential expression results.

## Function Description

The `analyze_pathway_volcano_custom` function generates a volcano plot highlighting genes from a specific pathway within differential expression results. It allows for flexible labeling and customization of the plot.

### Parameters

- `pathway_name`: Name of the pathway to analyze
- `gsea_results`: Results from Gene Set Enrichment Analysis (GSEA)
- `de_results`: Differential expression results
- `p_cutoff`: P-value cutoff for significance (default: 0.05)
- `fc_cutoff`: Fold change cutoff for significance (default: 2.0)
- `label_method`: Method for labeling genes ("default", "fc", "p", or "all")
- `max_overlaps`: Maximum number of label overlaps allowed (default: 100)

### Dependencies

The function requires the following R packages:

- ggplot2
- dplyr
- stringr
- ggrepel

## Usage

```R
source("analyze_pathway_volcano_custom.R")

# Assuming you have gsea_results and de_results
plot <- analyze_pathway_volcano_custom(
  pathway_name = "Your Pathway Name",
  gsea_results = gsea_results,
  de_results = de_results,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  label_method = "default",
  max_overlaps = 100
)
print(plot)
```

## Features

1. Extracts genes for the specified pathway from GSEA results
2. Prepares differential expression results with custom categories
3. Allows flexible gene labeling based on different criteria
4. Generates a customized volcano plot with:
   - Highlighted pathway genes
   - Significance thresholds for p-value and fold change
   - Custom color scheme
   - Repelled labels for selected genes

## Customization

You can customize the plot by adjusting the following parameters:

- `p_cutoff` and `fc_cutoff` for significance thresholds
- `label_method` for different gene labeling strategies:
  - "default": Label genes that are both in the pathway and significant
  - "fc": Label pathway genes based on fold change
  - "p": Label pathway genes based on p-value
  - "all": Label pathway genes that meet either fold change or p-value criteria
- `max_overlaps` to control label density

## Output

The function returns a ggplot object representing the volcano plot. The plot includes:

- Points representing genes, color-coded based on their significance and pathway membership
- Dashed lines indicating significance thresholds
- Labels for selected genes based on the chosen labeling method

## Contributing

Feel free to fork this repository and submit pull requests with improvements or bug fixes.

## License

[Specify your chosen license here]
