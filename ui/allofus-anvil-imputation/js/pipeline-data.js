/**
 * Pipeline-specific configuration for the 3 imputation products.
 *
 * To customize a pipeline, edit its entry below. Each pipeline has:
 *   - genomeOverviewHTML: left-side text in the reference panel section (HTML allowed)
 *   - totalGenomesCount / totalGenomesLabelHTML: donut chart callout
 *   - ancestryRows: table rows AND donut chart source of truth — count (formatted string),
 *       label, percent, and color (hex). The chart is derived from these; no separate array needed.
 *   - ancestryNoteHTML: footnote below the ancestry table (HTML allowed)
 *   - howItWorksSteps: ordered steps — title, bodyHTML (HTML allowed), img path, alt text
 *   - docsUrl: Documentation button link
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
        alt: "Install the command line tool",
      },
      {
        title: "Upload your data and launch",
        bodyHTML: "Upload your data to our secure environment and select parameters for your specific analysis.",
        img: "img/step3-data.png",
        alt: "Upload your data and launch",
      },
      {
        title: "Download your results",
        bodyHTML: "Download your results and review.",
        img: "img/step4-results.png",
        alt: "Download your results",
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
    name: "Low Pass WGS Imputation",
    tabDescription: "For low-coverage whole-genome sequencing data",
    pricePerSample: "$4.00",
    validationChart: {
      subtitle: "Imputation quality (r²) measured across 500 held-out samples sequenced at 0.5× coverage. Higher r² indicates more accurate genotype imputation.",
      xAxisLabel: "Minor Allele Frequency (MAF)",
      yAxisLabel: "Imputation Quality (r²)",
      labels: ["0.1–1%", "1–5%", "5–10%", "10–20%", "20–50%"],
      datasets: [
        {
          label: "All of Us + AnVIL Panel",
          data: [0.71, 0.89, 0.94, 0.97, 0.99],
          color: "#074770",
          dashed: false,
        },
        {
          label: "1000 Genomes",
          data: [0.45, 0.71, 0.84, 0.92, 0.97],
          color: "#ADB2BA",
          dashed: true,
        },
      ],
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
        alt: "Install the command line tool",
      },
      {
        title: "Upload your data and launch",
        bodyHTML: "Upload your data to our secure environment and select parameters for your specific analysis.",
        img: "img/step3-data.png",
        alt: "Upload your data and launch",
      },
      {
        title: "Download your results",
        bodyHTML: "Download your results and review.",
        img: "img/step4-results.png",
        alt: "Download your results",
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
    pricePerSample: "$1.00",
    genomeOverviewHTML: `The <i>All of Us</i> + AnVIL <br/>dataset contains <br/><span class="teal genome-count">515,000+ diverse <br/>genomes</span>`,
    totalGenomesCount: "13,000+",
    totalGenomesLabelHTML: `total genomes from <i>All of Us</i>`,
    ancestryRows: [
      { count:  "4,416", label: "European",              percent: "40%",  color: "#2A51B3" },
      { count:  "2,982", label: "African",               percent: "29%",  color: "#46A3E9" },
      { count:    "553", label: "Americas",              percent: "8%",   color: "#F6BD41" },
      { count:    "226", label: "East Asian",            percent: "3%",   color: "#5CC88D" },
      { count:  "1,710", label: "South Asian",           percent: "12%",  color: "#80C6EC" },
      { count:  "1,065", label: "Middle Eastern",        percent: "0.2%", color: "#ADB2BA" },
      { count:    "627", label: "Remaining Individuals", percent: "9%",   color: "#775FE5" },
    ],
    ancestryNoteHTML: `* Based on computed genetic ancestry on a combined dataset derived from the <i>All of Us</i> Curated Data Repository v8 release.`,
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
        alt: "Install the command line tool",
      },
      {
        title: "Upload your data and launch",
        bodyHTML: "Upload your data to our secure environment and select parameters for your specific analysis.",
        img: "img/step3-data.png",
        alt: "Upload your data and launch",
      },
      {
        title: "Download your results",
        bodyHTML: "Download your results and review.",
        img: "img/step4-results.png",
        alt: "Download your results",
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
};
