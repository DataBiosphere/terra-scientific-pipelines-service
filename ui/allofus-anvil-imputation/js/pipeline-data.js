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
    genomeOverviewHTML: `The <i>All of Us</i> + AnVIL <br/>dataset contains <br/><span class="teal genome-count">515,000+ diverse <br/>genomes</span>`,
    totalGenomesCount: "13,000+",
    totalGenomesLabelHTML: `total genomes from <i>All of Us</i> + AnVIL`,
    ancestryRows: [
      { count:  "4,416", label: "European",              percent: "40%",  color: "#2A51B3" },
      { count:  "2,982", label: "African",               percent: "29%",  color: "#46A3E9" },
      { count:    "553", label: "Americas",              percent: "8%",   color: "#F6BD41" },
      { count:    "226", label: "East Asian",            percent: "3%",   color: "#5CC88D" },
      { count:  "1,710", label: "South Asian",           percent: "12%",  color: "#80C6EC" },
      { count:  "1,065", label: "Middle Eastern",        percent: "0.2%", color: "#775FE5" },
      { count:    "627", label: "Remaining Individuals", percent: "9%",   color: "#ADB2BA" },
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
};
