/**
 * Pipeline-specific configuration for the 3 imputation products. Each pipeline's data lives in
 * its own file under js/pipelines/ and is collected here into PIPELINES, whose key order sets the
 * tab order. The objects aren't typed, but their schemas are documented below.
 */

/**
 * @typedef {Object} Pipeline
 * @property {string} name - required. Display name of the pipeline (e.g. "Array Imputation").
 * @property {string} tabDescription - required. Short description shown in the pipeline selection tab.
 * @property {number} priceForProfit - required. Cost per sample in dollars for for-profit users (e.g. 0.40).
 * @property {number} priceNonProfit - required. Discounted cost per sample in dollars for academic/nonprofit users (e.g. 0.30).
 * @property {number} maxNumSamples - required. Sample count at or above which the calculator directs users to contact us for alternative pricing instead of showing a Purchase button (e.g. 50000 means 49999 samples is the largest quantity priced by the calculator).
 * @property {string} genomeOverviewHTML - required. HTML string for the left-side text in the reference panel section.
 * @property {string} totalGenomesCount - required. Formatted count string for the donut chart callout (e.g. "515,000+").
 * @property {string} totalGenomesLabelHTML - required. HTML label beneath the donut chart count.
 * @property {AncestryRow[]} ancestryRows - required. Drives both the ancestry table rows and donut chart segments.
 * @property {string} ancestryNoteHTML - required. HTML footnote displayed below the ancestry table.
 * @property {HowItWorksStep[]} howItWorksSteps - required. Ordered steps shown in the "How It Works" section.
 * @property {string} docsUrl - required. URL for the Documentation button.
 * @property {string} pipelineKey - required. Pipeline identifier used for Mixpanel events, this page's `?pipeline=` query param, and the Teaspoon UI's `pipeline` query param (e.g. "array_imputation").
 * @property {ValidationChart[]} [validationCharts] - optional. Array of chart variants (e.g. SNP / INDEL); one toggle button is rendered per entry.
 * @property {ValidationPreprint} [validationPreprint] - optional. Preprint callout rendered below the validation chart.
 * @property {ComingSoon} [comingSoon] - optional. If present, the pipeline is shown as coming soon and all other fields are not required.
 */

/**
 * @typedef {Object} AncestryRow
 * @property {string} count - required. Formatted genome count (e.g. "254,416").
 * @property {string} label - required. Ancestry group label (e.g. "European").
 * @property {string} percent - required. Percentage of total (e.g. "49%").
 * @property {string} color - required. Hex color for the donut chart segment (e.g. "#2A51B3").
 */

/**
 * @typedef {Object} HowItWorksStep
 * @property {string} title - required. Step title.
 * @property {string} bodyHTML - required. HTML body text for the step.
 * @property {string} img - required. Path to the step illustration image.
 * @property {string} alt - required. Alt text for the step illustration image.
 */

/**
 * @typedef {Object} ValidationChart
 * @property {string} key - required. Unique identifier for the chart variant (e.g. "snp").
 * @property {string} buttonLabel - required. Label for the toggle button (e.g. "SNP").
 * @property {string} subtitle - required. Subtitle displayed above the chart.
 * @property {string} xAxisLabel - required. X-axis label.
 * @property {string} yAxisLabel - required. Y-axis label.
 * @property {string} [xAxisType] - optional. Chart.js axis type (e.g. "logarithmic"). Defaults to linear.
 * @property {Object.<number, string>} [tickLabels] - optional. Map of axis tick values to display labels, used when xAxisType is "logarithmic".
 * @property {ChartDataset[]} datasets - required for line charts. One or more datasets to plot on the chart.
 * @property {string} [chartType] - optional. "dumbbell" renders a grouped paired-marker chart (see the dumbbell fields below); omit for the default line chart.
 * @property {boolean} [animated] - optional. Set to false to skip Chart.js's draw-in animation when the chart first renders. Defaults to true.
 * @property {number} [yMin] - dumbbell only. Lower bound of the y-axis (e.g. 0.93).
 * @property {number} [yMax] - dumbbell only. Upper bound of the y-axis (e.g. 1.0).
 * @property {number} [yStepSize] - dumbbell only. Spacing between y-axis ticks (e.g. 0.01).
 * @property {number} [yTickDecimals] - dumbbell only. Decimal places for y-axis tick labels. Defaults to 2.
 * @property {string[]} [groups] - dumbbell only. Ordered x-axis group labels (e.g. superpopulations); each series must have one data point per group.
 * @property {string} [referenceLabel] - dumbbell only. Legend label for the open (reference) marker (e.g. "lrWGS panel").
 * @property {string} [imputedLabel] - dumbbell only. Legend label for the filled (imputed) marker (e.g. "srWGS imputed").
 * @property {string} [samplesLabel] - dumbbell only. Legend label for the per-sample dots. Defaults to "Individual sample"; shown only when any point has sample values.
 * @property {DumbbellSeries[]} [series] - dumbbell only. One entry per series fanned out within each group (e.g. SNV / INDEL / SV).
 */

/**
 * @typedef {Object} DumbbellSeries
 * @property {string} label - required. Legend label for the series (e.g. "SNV").
 * @property {string} color - required. Hex color for the series' legend swatch and range bands.
 * @property {DumbbellPoint[]} data - required. One point per entry in the chart's `groups`, in the same order.
 */

/**
 * @typedef {Object} DumbbellPoint
 * @property {number} reference - required. Value plotted as the open marker.
 * @property {number} imputed - required. Value plotted as the filled marker.
 * @property {number[]} [referenceSamples] - optional. Individual per-sample values for the reference side. When present (2+ values) they are drawn as small dots and as a half-violin (Gaussian KDE, Scott's bandwidth, truncated at the data extremes) to the left of the pair's centre line.
 * @property {number[]} [imputedSamples] - optional. Same as referenceSamples, for the imputed side (drawn to the right of the centre line).
 * @property {number[]} [referenceRange] - optional. [min, max] drawn as a shaded band behind the reference marker when no samples are given.
 * @property {number[]} [imputedRange] - optional. Same as referenceRange, for the imputed side.
 */

/**
 * @typedef {Object} ValidationPreprint
 * @property {string} badge - required. Text for the badge chip at the left of the pill (e.g. "New preprint").
 * @property {string} pillText - required. Main line of pill copy.
 * @property {string} linkText - required. Call-to-action text following the pill copy (e.g. "Read it on medRxiv").
 * @property {string} url - required. Preprint URL the pill links to.
 */

/**
 * @typedef {Object} ChartDataset
 * @property {string} label - required. Legend label for the dataset.
 * @property {{x: number, y: number}[]} data - required. Array of x/y data points.
 * @property {string} color - required. Hex color for the line (e.g. "#074770").
 * @property {boolean} dashed - required. Whether the line should be rendered as dashed.
 */

/**
 * @typedef {Object} ComingSoon
 * @property {string} message - required. Message to display on the coming soon banner.
 * @property {string} signupUrl - required. URL for the sign-up link.
 * @property {string} signupLabel - required. Display text for the sign-up link.
 */

const PIPELINES = {
  array: ARRAY_PIPELINE,
  lowpass: LOWPASS_PIPELINE,
  sv: SV_PIPELINE,
};
