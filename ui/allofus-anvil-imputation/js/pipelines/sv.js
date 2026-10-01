/**
 * SV Imputation pipeline data. See js/pipeline-data.js for the Pipeline schema and the PIPELINES registry.
 */
const SV_PIPELINE = {
  name: "SV Imputation",
  pipelineKey: "sv_imputation",
  tabDescription: "For structural variant inference",
  priceForProfit: 2.00,
  priceNonProfit: 1.50,
  maxNumSamples: 40000,
  validationCharts: [
    {
      key: "giab-easy",
      buttonLabel: "GIAB Easy Regions",
      chartType: "dumbbell",
      animated: false,
      subtitle: "F1 scores against assembly-based truth for SNVs, INDELs, and SVs on chromosome 20, in both the long-read panel and short-read held-out imputation, across 25 HPRC2/HGSVC3 samples stratified by superpopulation",
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
      buttonLabel: "GIAB Tandem Repeats & Homopolymers",
      chartType: "dumbbell",
      animated: false,
      subtitle: "F1 scores against assembly-based truth for SNVs, INDELs, and SVs on chromosome 20, in both the long-read panel and short-read held-out imputation, across 25 HPRC2/HGSVC3 samples stratified by superpopulation",
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
  ancestryNoteHTML: `* Based on computed genetic ancestry on a dataset derived from the <i>All of Us</i> Curated Data Repository v9 release.`,
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
};
