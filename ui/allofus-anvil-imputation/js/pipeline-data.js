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
    priceForProfit: 2.00,
    priceNonProfit: 1.50,
    // TODO: some values below are still placeholders (see "Placeholder" notes); replace with real SV pipeline data before launch
    maxNumSamples: 10000,
    // Computed from the vcfdist precision-recall summary (F1_SCORE, THRESHOLD=NONE) for 25 HPRC2/HGSVC3
    // samples on chr20, as provided in handoffs/plot2_data.json. Marker values are the per-group median;
    // the sample arrays hold each sample's F1 and drive the per-sample dots and the violins.
    validationCharts: [
      {
        key: "giab-easy",
        buttonLabel: "GIAB (Easy Regions)",
        chartType: "dumbbell",
        subtitle: "F1 score for SNVs, INDELs, and SVs on chromosome 20 in 25 HPRC2/HGSVC3 samples, comparing long-read reference panel calls with short-read imputed calls by superpopulation. Shaded bars show the per-sample range.",
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
              { // AFR (n=7)
                reference: 0.9918, imputed: 0.9890,
                referenceSamples: [0.991757, 0.991321, 0.990346, 0.992388, 0.99366, 0.989171, 0.992877],
                imputedSamples: [0.989505, 0.989889, 0.988805, 0.986686, 0.988967, 0.984329, 0.990455],
              },
              { // AMR (n=6)
                reference: 0.9906, imputed: 0.9901,
                referenceSamples: [0.993266, 0.993436, 0.989647, 0.990719, 0.980753, 0.990568],
                imputedSamples: [0.991935, 0.991279, 0.993528, 0.987843, 0.987852, 0.988907],
              },
              { // EAS (n=6)
                reference: 0.9911, imputed: 0.9834,
                referenceSamples: [0.994066, 0.988479, 0.992415, 0.99489, 0.989799, 0.985861],
                imputedSamples: [0.983045, 0.983766, 0.983958, 0.983467, 0.983248, 0.983422],
              },
              { // EUR (n=4)
                reference: 0.9904, imputed: 0.9889,
                referenceSamples: [0.983783, 0.994892, 0.990691, 0.990102],
                imputedSamples: [0.989733, 0.98765, 0.988766, 0.988956],
              },
              { // SAS (n=2)
                reference: 0.9926, imputed: 0.9843,
                referenceSamples: [0.994779, 0.990445],
                imputedSamples: [0.984252, 0.984399],
              },
            ],
          },
          {
            label: "INDEL",
            color: "#E0A322",
            data: [
              { // AFR (n=7)
                reference: 0.9653, imputed: 0.9622,
                referenceSamples: [0.96491, 0.965257, 0.963015, 0.969278, 0.968604, 0.964768, 0.967375],
                imputedSamples: [0.962211, 0.964195, 0.961335, 0.962465, 0.962137, 0.957541, 0.962181],
              },
              { // AMR (n=6)
                reference: 0.9625, imputed: 0.9638,
                referenceSamples: [0.968164, 0.967809, 0.96159, 0.963347, 0.960897, 0.960051],
                imputedSamples: [0.965949, 0.96505, 0.965385, 0.960663, 0.962452, 0.958723],
              },
              { // EAS (n=6)
                reference: 0.9625, imputed: 0.9530,
                referenceSamples: [0.965946, 0.942657, 0.965307, 0.966825, 0.959057, 0.959729],
                imputedSamples: [0.955459, 0.937003, 0.9559, 0.953249, 0.949606, 0.952758],
              },
              { // EUR (n=4)
                reference: 0.9616, imputed: 0.9588,
                referenceSamples: [0.957955, 0.966226, 0.961144, 0.962089],
                imputedSamples: [0.9568, 0.958128, 0.959402, 0.959515],
              },
              { // SAS (n=2)
                reference: 0.9606, imputed: 0.9538,
                referenceSamples: [0.962381, 0.95884],
                imputedSamples: [0.951408, 0.956132],
              },
            ],
          },
          {
            label: "SV",
            color: "#5CC88D",
            data: [
              { // AFR (n=7)
                reference: 0.9800, imputed: 0.9782,
                referenceSamples: [0.984127, 0.990566, 0.967818, 0.979976, 0.982488, 0.97308, 0.969422],
                imputedSamples: [0.978992, 0.976977, 0.963035, 0.979382, 0.978156, 0.97308, 0.978209],
              },
              { // AMR (n=6)
                reference: 0.9732, imputed: 0.9701,
                referenceSamples: [0.987952, 0.984615, 0.9774, 0.966197, 0.968962, 0.955286],
                imputedSamples: [0.970949, 0.974793, 0.971622, 0.954901, 0.969296, 0.948747],
              },
              { // EAS (n=6)
                reference: 0.9736, imputed: 0.9692,
                referenceSamples: [0.9875, 0.971979, 0.975274, 0.983889, 0.949489, 0.959805],
                imputedSamples: [0.969307, 0.966483, 0.969056, 0.973243, 0.949667, 0.979671],
              },
              { // EUR (n=4)
                reference: 0.9621, imputed: 0.9771,
                referenceSamples: [0.954149, 0.962442, 0.961849, 0.980392],
                imputedSamples: [0.988635, 0.967359, 0.961799, 0.986842],
              },
              { // SAS (n=2)
                reference: 0.9716, imputed: 0.9589,
                referenceSamples: [0.967612, 0.975595],
                imputedSamples: [0.962088, 0.955803],
              },
            ],
          },
        ],
      },
      {
        key: "giab-repeats",
        buttonLabel: "GIAB (Tandem Repeats & Homopolymers)",
        chartType: "dumbbell",
        subtitle: "F1 score for SNVs, INDELs, and SVs on chromosome 20 in 25 HPRC2/HGSVC3 samples, comparing long-read reference panel calls with short-read imputed calls by superpopulation. Shaded bars show the per-sample range.",
        xAxisLabel: "Superpopulation",
        yAxisLabel: "F1 Score",
        yMin: 0.74,
        yMax: 0.94,
        yStepSize: 0.02,
        groups: ["AFR", "AMR", "EAS", "EUR", "SAS"],
        referenceLabel: "Long-read WGS panel",
        imputedLabel: "Short-read WGS imputed",
        series: [
          {
            label: "SNV",
            color: "#2A51B3",
            data: [
              { // AFR (n=7)
                reference: 0.8120, imputed: 0.8037,
                referenceSamples: [0.78847, 0.812027, 0.80606, 0.812256, 0.832683, 0.798422, 0.830871],
                imputedSamples: [0.805898, 0.815433, 0.796085, 0.803732, 0.802606, 0.789591, 0.816922],
              },
              { // AMR (n=6)
                reference: 0.8065, imputed: 0.7997,
                referenceSamples: [0.825034, 0.82838, 0.773385, 0.762069, 0.79653, 0.816408],
                imputedSamples: [0.816089, 0.801309, 0.783006, 0.765668, 0.798294, 0.801069],
              },
              { // EAS (n=6)
                reference: 0.8067, imputed: 0.7757,
                referenceSamples: [0.829523, 0.775732, 0.814863, 0.831529, 0.798532, 0.785533],
                imputedSamples: [0.771088, 0.747871, 0.780403, 0.781154, 0.783411, 0.766352],
              },
              { // EUR (n=4)
                reference: 0.8088, imputed: 0.7882,
                referenceSamples: [0.819247, 0.818012, 0.769302, 0.799509],
                imputedSamples: [0.797221, 0.779255, 0.769716, 0.804471],
              },
              { // SAS (n=2)
                reference: 0.8080, imputed: 0.7656,
                referenceSamples: [0.816532, 0.799549],
                imputedSamples: [0.755532, 0.775634],
              },
            ],
          },
          {
            label: "INDEL",
            color: "#E0A322",
            data: [
              { // AFR (n=7)
                reference: 0.7901, imputed: 0.8105,
                referenceSamples: [0.78315, 0.790054, 0.780676, 0.8089, 0.842053, 0.774916, 0.812367],
                imputedSamples: [0.815751, 0.821939, 0.810514, 0.80299, 0.815855, 0.796562, 0.809872],
              },
              { // AMR (n=6)
                reference: 0.8000, imputed: 0.8275,
                referenceSamples: [0.795178, 0.830225, 0.788132, 0.802365, 0.797575, 0.804751],
                imputedSamples: [0.822348, 0.823376, 0.848307, 0.815554, 0.83769, 0.831668],
              },
              { // EAS (n=6)
                reference: 0.8111, imputed: 0.7894,
                referenceSamples: [0.86876, 0.791933, 0.798062, 0.864444, 0.817104, 0.80505],
                imputedSamples: [0.792006, 0.773248, 0.791733, 0.790702, 0.786906, 0.788143],
              },
              { // EUR (n=4)
                reference: 0.8004, imputed: 0.8115,
                referenceSamples: [0.803741, 0.864257, 0.796964, 0.795213],
                imputedSamples: [0.802781, 0.808761, 0.814258, 0.817382],
              },
              { // SAS (n=2)
                reference: 0.8309, imputed: 0.7861,
                referenceSamples: [0.869282, 0.79244],
                imputedSamples: [0.781064, 0.791227],
              },
            ],
          },
          {
            label: "SV",
            color: "#5CC88D",
            data: [
              { // AFR (n=7)
                reference: 0.8922, imputed: 0.8437,
                referenceSamples: [0.892176, 0.910098, 0.885631, 0.887142, 0.905899, 0.88816, 0.916027],
                imputedSamples: [0.845981, 0.853355, 0.824476, 0.840263, 0.843717, 0.831323, 0.854539],
              },
              { // AMR (n=6)
                reference: 0.9042, imputed: 0.8533,
                referenceSamples: [0.906158, 0.911286, 0.891075, 0.902311, 0.868717, 0.923285],
                imputedSamples: [0.887851, 0.859288, 0.868414, 0.840341, 0.816059, 0.84735],
              },
              { // EAS (n=6)
                reference: 0.8843, imputed: 0.8055,
                referenceSamples: [0.892242, 0.848319, 0.876291, 0.892875, 0.901181, 0.866152],
                imputedSamples: [0.812964, 0.763504, 0.81983, 0.798045, 0.821673, 0.797431],
              },
              { // EUR (n=4)
                reference: 0.9005, imputed: 0.8342,
                referenceSamples: [0.901465, 0.884271, 0.902988, 0.899455],
                imputedSamples: [0.84014, 0.836229, 0.827995, 0.832109],
              },
              { // SAS (n=2)
                reference: 0.8775, imputed: 0.7929,
                referenceSamples: [0.878175, 0.876789],
                imputedSamples: [0.78402, 0.801774],
              },
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
    ancestryNoteHTML: `* Based on computed genetic ancestry on a dataset derived from the <i>All of Us</i> Curated Data Repository v8 release.`,
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
