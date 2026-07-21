/**
 * Pipeline-specific configuration for the 3 imputation products. This isn't typed but
 * the object schemas are shown below to help construct the necessary pipeline data.
 */

/**
 * @typedef {Object} Pipeline
 * @property {string} name - required. Display name of the pipeline (e.g. "Array Imputation").
 * @property {string} tabDescription - required. Short description shown in the pipeline selection tab.
 * @property {string} pricePerSample - required. Cost per sample (e.g. "$0.40").
 * @property {string} genomeOverviewHTML - required. HTML string for the left-side text in the reference panel section.
 * @property {string} totalGenomesCount - required. Formatted count string for the donut chart callout (e.g. "515,000+").
 * @property {string} totalGenomesLabelHTML - required. HTML label beneath the donut chart count.
 * @property {AncestryRow[]} ancestryRows - required. Drives both the ancestry table rows and donut chart segments.
 * @property {string} ancestryNoteHTML - required. HTML footnote displayed below the ancestry table.
 * @property {HowItWorksStep[]} howItWorksSteps - required. Ordered steps shown in the "How It Works" section.
 * @property {string} docsUrl - required. URL for the Documentation button.
 * @property {ValidationChart[]} [validationCharts] - optional. Array of chart variants (e.g. SNP / INDEL); one toggle button is rendered per entry.
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
 * @property {ChartDataset[]} datasets - required. One or more datasets to plot on the chart.
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
  array: {
    name: "Array Imputation",
    tabDescription: "For array-based genotype data",
    pricePerSample: "$0.40",
    genomeOverviewHTML: `The <i>All of Us</i> + AnVIL <br/>dataset contains <br/><span class="teal genome-count">515,000+ diverse <br/>genomes</span>`,
    totalGenomesCount: "515,000+",
    totalGenomesLabelHTML: `total genomes from <i>All of Us</i> + AnVIL`,
    ancestryRows: [
      { count: "254,416", label: "European",              percent: "49%",  color: "#2A51B3" },
      { count: "101,982", label: "African",               percent: "20%",  color: "#46A3E9" },
      { count:  "90,553", label: "Americas",              percent: "18%",  color: "#F6BD41" },
      { count:  "13,226", label: "East Asian",            percent: "3%",   color: "#80C6EC" },
      { count:   "9,710", label: "South Asian",           percent: "2%",   color: "#775FE5" },
      { count:   "1,065", label: "Middle Eastern",        percent: "0.2%", color: "#ADB2BA" },
      { count:  "44,627", label: "Remaining Individuals", percent: "9%",   color: "#5CC88D" },
    ],
    ancestryNoteHTML: `* Based on computed genetic ancestry on a combined dataset derived from the <i>All of Us</i> Curated Data Repository v8 release and AnVIL Centers for Common Disease Genomics.`,
    howItWorksSteps: [
      {
        title: "Create an account",
        bodyHTML: "Create a Terra account to get started.",
        img: "img/step1-create-account.png",
        alt: "Create an account",
      },
      {
        title: "Pick your preferred method",
        bodyHTML: `Visit our <a href="https://services.terra.bio/" target="_blank">web interface</a> or <a href="https://broadscientificservices.zendesk.com/hc/en-us/articles/39901313672859" target="_blank">install our command-line tool</a> in your preferred environment.`,
        img: "img/step2-download.png",
        alt: "Install the command line tool or use the web-based UI",
      },
      {
        title: "Bring your data and launch",
        bodyHTML: "Upload your data to our secure environment or use existing data in Google Cloud Storage and select parameters for your specific analysis.",
        img: "img/step3-data.png",
        alt: "Bring your data and launch",
      },
      {
        title: "Retrieve your results",
        bodyHTML: "Download your results or have them delivered to a Google Cloud Storage bucket of your choice.",
        img: "img/step4-results.png",
        alt: "Retrieve your results",
      },
      {
        title: "Run analyses on the imputed dataset",
        bodyHTML: `Optionally continue your analysis in Terra where we provide genome-wide analysis and polygenic risk score pipelines to enable downstream analysis on the imputed dataset.`,
        img: "img/step5-analysis.png",
        alt: "Run analyses on the imputed dataset",
      },
    ],
    docsUrl: "https://broadscientificservices.zendesk.com/hc/en-us/categories/39900993442459",
  },

  lowpass: {
    name: "Low-Pass WGS Imputation",
    tabDescription: "For low-coverage whole-genome sequencing data",
    pricePerSample: "$4.00",
    validationCharts: [
      {
        key: "snp",
        buttonLabel: "SNP",
        subtitle: "Aggregate r² for SNPs and INDELs across 42 held-out samples, benchmarked against 1000 Genomes",
        xAxisLabel: "Allele Frequency (AF)",
        yAxisLabel: "Imputation Quality (r²)",
        xAxisType: "logarithmic",
        tickLabels: { 0.001: "0.001", 0.01: "0.01", 0.1: "0.1", 1: "1" },
        datasets: [
          {
            label: "All of Us + AnVIL",
            data: [
              { x: 0.00025,   y: 0.7979657967470237 },
              { x: 0.00067,   y: 0.9075451574147143 },
              { x: 0.001202,  y: 0.9263276644021905 },
              { x: 0.002157,  y: 0.942691515927857 },
              { x: 0.00387,   y: 0.9576006742544999 },
              { x: 0.006944,  y: 0.9693611579737381 },
              { x: 0.012461,  y: 0.9791937725752142 },
              { x: 0.022361,  y: 0.9858736919972143 },
              { x: 0.040125,  y: 0.9907515842588095 },
              { x: 0.072001,  y: 0.9937542064404286 },
              { x: 0.1292,    y: 0.9953788671581429 },
              { x: 0.231839,  y: 0.995972844717881 },
              { x: 0.416018,  y: 0.9959094803146668 },
              { x: 0.746513,  y: 0.9964782864694047 },
            ],
            color: "#074770",
            dashed: false,
          },
          {
            label: "1000 Genomes",
            data: [
              { x: 0.00025,   y: 0.46706493557885714 },
              { x: 0.00067,   y: 0.7092035360227142 },
              { x: 0.001202,  y: 0.7771167021560953 },
              { x: 0.002157,  y: 0.837234889633381 },
              { x: 0.00387,   y: 0.8879356314405953 },
              { x: 0.006944,  y: 0.9241190264634525 },
              { x: 0.012461,  y: 0.9499619409738096 },
              { x: 0.022361,  y: 0.9672541219321905 },
              { x: 0.040125,  y: 0.9791080135601428 },
              { x: 0.072001,  y: 0.9861876462721666 },
              { x: 0.1292,    y: 0.9894929043439524 },
              { x: 0.231839,  y: 0.9907366901243572 },
              { x: 0.416018,  y: 0.991374024785643 },
              { x: 0.746513,  y: 0.992128286392619 },
            ],
            color: "#ADB2BA",
            dashed: true,
          },
        ],
      },
      {
        key: "indel",
        buttonLabel: "INDEL",
        subtitle: "Aggregate r² for SNPs and INDELs across 42 held-out samples, benchmarked against 1000 Genomes",
        xAxisLabel: "Allele Frequency (AF)",
        yAxisLabel: "Imputation Quality (r²)",
        xAxisType: "logarithmic",
        tickLabels: { 0.001: "0.001", 0.01: "0.01", 0.1: "0.1", 1: "1" },
        datasets: [
          {
            label: "All of Us + AnVIL",
            data: [
              { x: 0.00025,   y: 0.5883368902868095 },
              { x: 0.00067,   y: 0.6801708592414761 },
              { x: 0.001202,  y: 0.7175159369097619 },
              { x: 0.002157,  y: 0.7562319707075 },
              { x: 0.00387,   y: 0.7959195920611667 },
              { x: 0.006944,  y: 0.8295353032107619 },
              { x: 0.012461,  y: 0.8641983264370475 },
              { x: 0.022361,  y: 0.8944237852621191 },
              { x: 0.040125,  y: 0.9215695944036428 },
              { x: 0.072001,  y: 0.9452490735752144 },
              { x: 0.1292,    y: 0.9618026830943333 },
              { x: 0.231839,  y: 0.972471692681119 },
              { x: 0.416018,  y: 0.980283187208619 },
              { x: 0.746513,  y: 0.9846723969148333 },
            ],
            color: "#074770",
            dashed: false,
          },
          {
            label: "1000 Genomes",
            data: [
              { x: 0.00025,   y: 0.2303683394470238 },
              { x: 0.00067,   y: 0.3252659980248572 },
              { x: 0.001202,  y: 0.3884575066018809 },
              { x: 0.002157,  y: 0.46786471277080954 },
              { x: 0.00387,   y: 0.5518749506321191 },
              { x: 0.006944,  y: 0.6174154834853334 },
              { x: 0.012461,  y: 0.6874531574870952 },
              { x: 0.022361,  y: 0.7426582641191904 },
              { x: 0.040125,  y: 0.795741420979619 },
              { x: 0.072001,  y: 0.8438215520944762 },
              { x: 0.1292,    y: 0.8870549570013571 },
              { x: 0.231839,  y: 0.921170141711262 },
              { x: 0.416018,  y: 0.9487033790830476 },
              { x: 0.746513,  y: 0.9642887887423333 },
            ],
            color: "#ADB2BA",
            dashed: true,
          },
        ],
      },
    ],
    genomeOverviewHTML: `The <i>All of Us</i> + AnVIL <br/>dataset contains <br/><span class="teal genome-count">515,000+ diverse <br/>genomes</span>`,
    totalGenomesCount: "515,000+",
    totalGenomesLabelHTML: `total genomes from <i>All of Us</i> + AnVIL`,
    ancestryRows: [
      { count: "254,416", label: "European",              percent: "49%",  color: "#2A51B3" },
      { count: "101,982", label: "African",               percent: "20%",  color: "#46A3E9" },
      { count:  "90,553", label: "Americas",              percent: "18%",  color: "#F6BD41" },
      { count:  "13,226", label: "East Asian",            percent: "3%",   color: "#80C6EC" },
      { count:   "9,710", label: "South Asian",           percent: "2%",   color: "#775FE5" },
      { count:   "1,065", label: "Middle Eastern",        percent: "0.2%", color: "#ADB2BA" },
      { count:  "44,627", label: "Remaining Individuals", percent: "9%",   color: "#5CC88D" },
    ],
    ancestryNoteHTML: `* Based on computed genetic ancestry on a combined dataset derived from the <i>All of Us</i> Curated Data Repository v8 release and AnVIL Centers for Common Disease Genomics.`,
    howItWorksSteps: [
      {
        title: "Create an account",
        bodyHTML: "Create a Terra account to get started.",
        img: "img/step1-create-account.png",
        alt: "Create an account",
      },
      {
        title: "Pick your preferred method",
        bodyHTML: `Visit our <a href="https://services.terra.bio/" target="_blank">web interface</a> or <a href="https://broadscientificservices.zendesk.com/hc/en-us/articles/39901313672859" target="_blank">install our command-line tool</a> in your preferred environment.`,
        img: "img/step2-download.png",
        alt: "Install the command line tool or use the web-based UI",
      },
      {
        title: "Bring your data and launch",
        bodyHTML: "Bring your Google Cloud-hosted data to our secure environment and select parameters for your specific analysis.",
        img: "img/step3-data.png",
        alt: "Bring your data and launch",
      },
      {
        title: "Retrieve your results",
        bodyHTML: "Download your results or have them delivered to a Google Cloud Storage bucket of your choice.",
        img: "img/step4-results.png",
        alt: "Retrieve your results",
      },
      {
        title: "Run analyses on the imputed dataset",
        bodyHTML: `Optionally continue your analysis in Terra where we provide genome-wide analysis and polygenic risk score pipelines to enable downstream analysis on the imputed dataset.`,
        img: "img/step5-analysis.png",
        alt: "Run analyses on the imputed dataset",
      },
    ],
    docsUrl: "https://broadscientificservices.zendesk.com/hc/en-us/categories/39900993442459",
  },

  sv: {
    name: "SV Imputation",
    tabDescription: "For structural variant inference",
    comingSoon: {
      message: "We're launching a new SV Imputation pipeline in late 2026.",
      signupUrl: "#",
      signupLabel: "Sign up to be notified",
    },
  },
};
