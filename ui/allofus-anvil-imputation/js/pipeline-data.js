/**
 * Pipeline-specific configuration for the 3 imputation products. This isn't typed but
 * the object schemas are shown below to help construct the necessary pipeline data.
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
 * @property {string} pipelineKey - required. Pipeline identifier used for Mixpanel events and the Teaspoon UI's `pipeline` query param (e.g. "array_imputation").
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
 * @property {number} [yMin] - dumbbell only. Lower bound of the y-axis (e.g. 0.93).
 * @property {number} [yMax] - dumbbell only. Upper bound of the y-axis (e.g. 1.0).
 * @property {number} [yStepSize] - dumbbell only. Spacing between y-axis ticks (e.g. 0.01).
 * @property {number} [yTickDecimals] - dumbbell only. Decimal places for y-axis tick labels. Defaults to 2.
 * @property {string[]} [groups] - dumbbell only. Ordered x-axis group labels (e.g. superpopulations); each series must have one data point per group.
 * @property {string} [referenceLabel] - dumbbell only. Legend label for the open (reference) marker (e.g. "lrWGS panel").
 * @property {string} [imputedLabel] - dumbbell only. Legend label for the filled (imputed) marker (e.g. "srWGS imputed").
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
 * @property {number[]} [referenceRange] - optional. [min, max] drawn as a shaded band behind the reference marker.
 * @property {number[]} [imputedRange] - optional. [min, max] drawn as a shaded band behind the imputed marker.
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
  array: {
    name: "Array Imputation",
    pipelineKey: "array_imputation",
    tabDescription: "For array-based genotype data",
    priceForProfit: 0.40,
    priceNonProfit: 0.30,
    maxNumSamples: 50000,
    validationCharts: [
      {
        key: "snp",
        buttonLabel: "SNP",
        subtitle: "Aggregate r² for SNPs and INDELs across 42 held-out samples, from diverse reference populations, benchmarked against TOPMed",
        xAxisLabel: "Allele Frequency (AF)",
        yAxisLabel: "Imputation Quality (r²)",
        xAxisType: "logarithmic",
        tickLabels: { 0.001: "0.001", 0.01: "0.01", 0.1: "0.1", 1: "1" },
        datasets: [
          {
            label: "All of Us + AnVIL",
            data: [
              { x: 0.00025,   y: 0.6501123264783095 },
              { x: 0.00067,   y: 0.7804926872044524 },
              { x: 0.001202,  y: 0.812926018418381 },
              { x: 0.002157,  y: 0.8456072153253096 },
              { x: 0.00387,   y: 0.8768931761803332 },
              { x: 0.006944,  y: 0.9027232215536666 },
              { x: 0.012461,  y: 0.9311253899843334 },
              { x: 0.022361,  y: 0.9500231940289287 },
              { x: 0.040125,  y: 0.9662221480160239 },
              { x: 0.072001,  y: 0.9774507448003572 },
              { x: 0.1292,    y: 0.9838099675994763 },
              { x: 0.231839,  y: 0.9867967533926906 },
              { x: 0.416018,  y: 0.9883669655030238 },
              { x: 0.746513,  y: 0.9894291085064763 },
            ],
            color: "#074770",
            dashed: false,
          },
          {
            label: "TOPMed",
            data: [
              { x: 0.00025,   y: 0.5061692868880951 },
              { x: 0.00067,   y: 0.6771790604641905 },
              { x: 0.001202,  y: 0.7247445797786904 },
              { x: 0.002157,  y: 0.7716195364405716 },
              { x: 0.00387,   y: 0.818560953817619 },
              { x: 0.006944,  y: 0.8581645502408334 },
              { x: 0.012461,  y: 0.8969946711791428 },
              { x: 0.022361,  y: 0.9255383077945476 },
              { x: 0.040125,  y: 0.9517807583236666 },
              { x: 0.072001,  y: 0.967772856609119 },
              { x: 0.1292,    y: 0.9768377178131905 },
              { x: 0.231839,  y: 0.9807297596859763 },
              { x: 0.416018,  y: 0.9831181136565238 },
              { x: 0.746513,  y: 0.9847500581101666 },
            ],
            color: "#ADB2BA",
            dashed: true,
          },
        ],
      },
      {
        key: "indel",
        buttonLabel: "INDEL",
        subtitle: "Aggregate r² for SNPs and INDELs across 42 held-out samples, from diverse reference populations, benchmarked against TOPMed",
        xAxisLabel: "Allele Frequency (AF)",
        yAxisLabel: "Imputation Quality (r²)",
        xAxisType: "logarithmic",
        tickLabels: { 0.001: "0.001", 0.01: "0.01", 0.1: "0.1", 1: "1" },
        datasets: [
          {
            label: "All of Us + AnVIL",
            data: [
              { x: 0.00025,   y: 0.6108649517158333 },
              { x: 0.00067,   y: 0.7402963440531904 },
              { x: 0.001202,  y: 0.7830736015996429 },
              { x: 0.002157,  y: 0.8213534964221191 },
              { x: 0.00387,   y: 0.8623922362660714 },
              { x: 0.006944,  y: 0.8951441141059286 },
              { x: 0.012461,  y: 0.9270198033703094 },
              { x: 0.022361,  y: 0.9449179827933809 },
              { x: 0.040125,  y: 0.9623308255877618 },
              { x: 0.072001,  y: 0.975082363039381 },
              { x: 0.1292,    y: 0.9818741673537618 },
              { x: 0.231839,  y: 0.985803179004262 },
              { x: 0.416018,  y: 0.9876879168465 },
              { x: 0.746513,  y: 0.9890449802521429 },
            ],
            color: "#074770",
            dashed: false,
          },
          {
            label: "TOPMed",
            data: [
              { x: 0.00025,   y: 0.47002659892469045 },
              { x: 0.00067,   y: 0.6422049021947381 },
              { x: 0.001202,  y: 0.6919790244210239 },
              { x: 0.002157,  y: 0.7485425035227857 },
              { x: 0.00387,   y: 0.8052055125385952 },
              { x: 0.006944,  y: 0.8448814384093571 },
              { x: 0.012461,  y: 0.8890615502314048 },
              { x: 0.022361,  y: 0.915717333963119 },
              { x: 0.040125,  y: 0.9422920466460953 },
              { x: 0.072001,  y: 0.9591011970989285 },
              { x: 0.1292,    y: 0.9678963193053095 },
              { x: 0.231839,  y: 0.9708382792160952 },
              { x: 0.416018,  y: 0.969575388000881 },
              { x: 0.746513,  y: 0.9433954912202619 },
            ],
            color: "#ADB2BA",
            dashed: true,
          },
        ],
      },
    ],
    validationPreprint: {
      badge: "New preprint",
      pillText: "The reference panel and these benchmarks are described in our preprint",
      linkText: "Read it on medRxiv",
      url: "https://www.medrxiv.org/content/10.64898/2026.08.25.26361247v1",
    },
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
        bodyHTML: `Download your results or have them delivered to a Google Cloud Storage bucket of your choice.<div class="how-step-note"><svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"/><polyline points="12 6 12 12 16 14"/></svg> Amazon S3 support is coming soon.</div>`,
        img: "img/step4-results.png",
        alt: "Retrieve your results",
      }
    ],
    docsUrl: "https://broadscientificservices.zendesk.com",
  },

  lowpass: {
    name: "Low-Pass WGS Imputation",
    pipelineKey: "low_pass_imputation",
    tabDescription: "For low-coverage whole-genome sequencing data",
    priceForProfit: 4.00,
    priceNonProfit: 3.50,
    maxNumSamples: 10000,
    validationCharts: [
      {
        key: "snp",
        buttonLabel: "SNP",
        subtitle: "Aggregate r² for SNPs and INDELs across 42 held-out samples, from diverse reference populations, benchmarked against 1000 Genomes",
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
        subtitle: "Aggregate r² for SNPs and INDELs across 42 held-out samples, from diverse reference populations, benchmarked against 1000 Genomes",
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
        bodyHTML: `Download your results or have them delivered to a Google Cloud Storage bucket of your choice.<div class="how-step-note"><svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"/><polyline points="12 6 12 12 16 14"/></svg> Amazon S3 support is coming soon.</div>`,
        img: "img/step4-results.png",
        alt: "Retrieve your results",
      }
    ],
    docsUrl: "https://broadscientificservices.zendesk.com",
  },

  sv: {
    name: "SV Imputation",
    pipelineKey: "sv_imputation",
    tabDescription: "For structural variant inference",
    // TODO: placeholder values below, replace with real SV pipeline data before launch
    priceForProfit: 0.00,
    priceNonProfit: 0.00,
    maxNumSamples: 10000,
    // Values were read from the handoff figure (handoffs/25_HPRC_HGSVC_chr20_F1.png). Marker values are
    // the plotted means; ranges are the vertical extent of each half-violin (per-sample min/max).
    // TODO: replace with the source numbers behind that figure when available.
    validationCharts: [
      {
        key: "giab-easy",
        buttonLabel: "GIAB easy regions",
        chartType: "dumbbell",
        subtitle: "F1 score for SNVs, INDELs, and SVs on chr20 in 25 HPRC2/HGSVC3 samples, comparing long-read (lrWGS) reference panel calls with short-read (srWGS) imputed calls by superpopulation. Shaded bars show the per-sample range.",
        xAxisLabel: "Superpopulation",
        yAxisLabel: "F1 Score",
        yMin: 0.93,
        yMax: 1.0,
        yStepSize: 0.01,
        groups: ["AFR", "AMR", "EAS", "EUR", "SAS"],
        referenceLabel: "Long-read WGS panel",
        imputedLabel: "Short-read WGS imputed",
        series: [
          {
            label: "SNV",
            color: "#2A51B3",
            data: [
              { reference: 0.9917, imputed: 0.9889, referenceRange: [0.9892, 0.9936], imputedRange: [0.9844, 0.9904] }, // AFR
              { reference: 0.9906, imputed: 0.9901, referenceRange: [0.9808, 0.9934], imputedRange: [0.9879, 0.9935] }, // AMR
              { reference: 0.9911, imputed: 0.9834, referenceRange: [0.9859, 0.9948], imputedRange: [0.9833, 0.9834] }, // EAS
              { reference: 0.9904, imputed: 0.9888, referenceRange: [0.9838, 0.9948], imputedRange: [0.9877, 0.9897] }, // EUR
              { reference: 0.9926, imputed: 0.9843, referenceRange: [0.9905, 0.9947], imputedRange: [0.9843, 0.9843] }, // SAS
            ],
          },
          {
            label: "INDEL",
            color: "#E0A322",
            data: [
              { reference: 0.9652, imputed: 0.9622, referenceRange: [0.9630, 0.9692], imputedRange: [0.9576, 0.9641] }, // AFR
              { reference: 0.9624, imputed: 0.9637, referenceRange: [0.9601, 0.9681], imputedRange: [0.9587, 0.9659] }, // AMR
              { reference: 0.9625, imputed: 0.9530, referenceRange: [0.9426, 0.9668], imputedRange: [0.9371, 0.9559] }, // EAS
              { reference: 0.9616, imputed: 0.9587, referenceRange: [0.9580, 0.9662], imputedRange: [0.9568, 0.9594] }, // EUR
              { reference: 0.9606, imputed: 0.9538, referenceRange: [0.9588, 0.9624], imputedRange: [0.9514, 0.9561] }, // SAS
            ],
          },
          {
            label: "SV",
            color: "#5CC88D",
            data: [
              { reference: 0.9800, imputed: 0.9781, referenceRange: [0.9678, 0.9905], imputedRange: [0.9630, 0.9793] }, // AFR
              { reference: 0.9731, imputed: 0.9701, referenceRange: [0.9553, 0.9879], imputedRange: [0.9488, 0.9748] }, // AMR
              { reference: 0.9736, imputed: 0.9692, referenceRange: [0.9495, 0.9875], imputedRange: [0.9498, 0.9796] }, // EAS
              { reference: 0.9621, imputed: 0.9771, referenceRange: [0.9542, 0.9803], imputedRange: [0.9618, 0.9886] }, // EUR
              { reference: 0.9716, imputed: 0.9589, referenceRange: [0.9676, 0.9755], imputedRange: [0.9558, 0.9620] }, // SAS
            ],
          },
        ],
      },
      {
        key: "giab-repeats",
        buttonLabel: "GIAB tandem repeats & homopolymers",
        chartType: "dumbbell",
        subtitle: "F1 score for SNVs, INDELs, and SVs on chr20 in 25 HPRC2/HGSVC3 samples, comparing long-read (lrWGS) reference panel calls with short-read (srWGS) imputed calls by superpopulation. Shaded bars show the per-sample range.",
        xAxisLabel: "Superpopulation",
        yAxisLabel: "F1 Score",
        yMin: 0.74,
        yMax: 0.94,
        yStepSize: 0.02,
        groups: ["AFR", "AMR", "EAS", "EUR", "SAS"],
        referenceLabel: "lrWGS panel",
        imputedLabel: "srWGS imputed",
        series: [
          {
            label: "SNV",
            color: "#2A51B3",
            data: [
              { reference: 0.8120, imputed: 0.8037, referenceRange: [0.7886, 0.8325], imputedRange: [0.7897, 0.8167] }, // AFR
              { reference: 0.8065, imputed: 0.7997, referenceRange: [0.7622, 0.8282], imputedRange: [0.7658, 0.8159] }, // AMR
              { reference: 0.8068, imputed: 0.7758, referenceRange: [0.7759, 0.8314], imputedRange: [0.7481, 0.7834] }, // EAS
              { reference: 0.8087, imputed: 0.7882, referenceRange: [0.7695, 0.8191], imputedRange: [0.7699, 0.8044] }, // EUR
              { reference: 0.8080, imputed: 0.7657, referenceRange: [0.7997, 0.8164], imputedRange: [0.7557, 0.7755] }, // SAS
            ],
          },
          {
            label: "INDEL",
            color: "#E0A322",
            data: [
              { reference: 0.7900, imputed: 0.8105, referenceRange: [0.7751, 0.8419], imputedRange: [0.7967, 0.8218] }, // AFR
              { reference: 0.8000, imputed: 0.8275, referenceRange: [0.7882, 0.8297], imputedRange: [0.8156, 0.8481] }, // AMR
              { reference: 0.8110, imputed: 0.7895, referenceRange: [0.7921, 0.8686], imputedRange: [0.7737, 0.7919] }, // EAS
              { reference: 0.8004, imputed: 0.8115, referenceRange: [0.7953, 0.8641], imputedRange: [0.8029, 0.8173] }, // EUR
              { reference: 0.8308, imputed: 0.7861, referenceRange: [0.7925, 0.8692], imputedRange: [0.7812, 0.7911] }, // SAS
            ],
          },
          {
            label: "SV",
            color: "#5CC88D",
            data: [
              { reference: 0.8920, imputed: 0.8437, referenceRange: [0.8856, 0.9158], imputedRange: [0.8246, 0.8544] }, // AFR
              { reference: 0.9040, imputed: 0.8533, referenceRange: [0.8688, 0.9230], imputedRange: [0.8162, 0.8877] }, // AMR
              { reference: 0.8841, imputed: 0.8055, referenceRange: [0.8484, 0.9010], imputedRange: [0.7637, 0.8216] }, // EAS
              { reference: 0.9003, imputed: 0.8342, referenceRange: [0.8843, 0.9028], imputedRange: [0.8281, 0.8400] }, // EUR
              { reference: 0.8773, imputed: 0.7929, referenceRange: [0.8773, 0.8775], imputedRange: [0.7842, 0.8016] }, // SAS
            ],
          },
        ],
      },
    ],
    genomeOverviewHTML: `The <i>All of Us</i><br/>dataset contains <br/><span class="teal genome-count">12,500+ diverse <br/>genomes</span>`,
    totalGenomesCount: "12,500+",
    totalGenomesLabelHTML: `total genomes from <i>All of Us</i>`,
    ancestryRows: [
      { count: "2,744", label: "African",               percent: "22%", color: "#46A3E9" },
      { count: "2,131", label: "Americas",              percent: "17%", color: "#2A51B3" },
      { count: "1,901", label: "European",              percent: "15%", color: "#F6BD41" },
      { count: "1,373", label: "East Asian",            percent: "11%", color: "#80C6EC" },
      { count: "1,241", label: "South Asian",           percent: "10%", color: "#775FE5" },
      { count: "416",   label: "Middle Eastern",        percent: "3%",  color: "#ADB2BA" },
      { count: "2,748", label: "Remaining Individuals", percent: "22%", color: "#5CC88D" },
    ],
    ancestryNoteHTML: `* Based on computed genetic ancestry on a combined dataset derived from the All of Us Curated Data Repository v8 release.`,
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
        bodyHTML: "Placeholder: describe SV input data requirements and parameters.",
        img: "img/step3-data.png",
        alt: "Bring your data and launch",
      },
      {
        title: "Retrieve your results",
        bodyHTML: `Download your results or have them delivered to a Google Cloud Storage bucket of your choice.<div class="how-step-note"><svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"/><polyline points="12 6 12 12 16 14"/></svg> Amazon S3 support is coming soon.</div>`,
        img: "img/step4-results.png",
        alt: "Retrieve your results",
      }
    ],
    docsUrl: "https://broadscientificservices.zendesk.com",
  },
};
