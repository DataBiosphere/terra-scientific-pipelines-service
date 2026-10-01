/**
 * Scientific validation section (pipeline-specific, optional — hidden when p.validationCharts is absent).
 * Renders one toggle button per entry in p.validationCharts (e.g. SNP / INDEL), plus an optional
 * preprint pill below the chart (p.validationPreprint).
 *
 * The chart itself is drawn by the type named in vc.chartType (default "line"), looked up in
 * VALIDATION_CHART_TYPES — see js/components/chart-common.js and the js/components/*-chart.js files.
 */
let _validationChart = null;

const VALIDATION_ARROW_SVG = `<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round" style="flex-shrink:0" aria-hidden="true">
  <line x1="5" y1="12" x2="19" y2="12"/><polyline points="12 5 19 12 12 19"/></svg>`;

function validationPreprintHTML(pp) {
  if (!pp) return '';
  return `
    <a class="validation-preprint" href="${pp.url}" target="_blank" rel="noopener">
      <div class="validation-preprint-copy">
        <span class="validation-preprint-badge">${pp.badge}</span>
        <span class="validation-preprint-text">${pp.pillText}</span>
      </div>
      <span class="validation-preprint-link">${pp.linkText} ${VALIDATION_ARROW_SVG}</span>
    </a>`;
}

function renderValidationChart(vc) {
  const canvas = document.getElementById('validationChartCanvas');
  if (_validationChart) { _validationChart.destroy(); _validationChart = null; }

  Chart.defaults.font.family = 'Montserrat';

  _validationChart = new Chart(canvas.getContext('2d'), validationChartType(vc).config(vc));
}

function validationChartType(vc) {
  const type = VALIDATION_CHART_TYPES[vc.chartType || 'line'];
  if (!type) throw new Error(`Unknown validation chart type: ${vc.chartType}`);
  return type;
}

function renderValidationSection(p) {
  const container = document.getElementById('frame-validation');
  const charts = p.validationCharts;
  if (!charts || !charts.length) {
    container.innerHTML = '';
    container.style.display = 'none';
    if (_validationChart) { _validationChart.destroy(); _validationChart = null; }
    return;
  }
  container.style.display = '';

  let activeKey = charts[0].key;

  const renderChartToggle = () => {
    if (charts.length < 2) return '';
    return `<div class="validation-chart-toggle">
      ${charts.map(c => `<button class="chart-toggle-btn${c.key === activeKey ? ' active' : ''}" data-chart-key="${c.key}">${c.buttonLabel}</button>`).join('')}
    </div>`;
  };

  const draw = () => {
    const vc = charts.find(c => c.key === activeKey);
    container.innerHTML = `
      <div class="validation-header">
        Has the imputation service been scientifically validated?
        <div class="validation-subtext">${vc.subtitle}</div>
      </div>
      ${renderChartToggle()}
      <div class="validation-chart-wrapper">
        ${validationChartType(vc).legendHTML(vc)}
        <canvas id="validationChartCanvas"></canvas>
      </div>
      ${validationPreprintHTML(p.validationPreprint)}`;

    container.querySelectorAll('.chart-toggle-btn').forEach(btn => {
      btn.addEventListener('click', () => {
        if (btn.dataset.chartKey === activeKey) return;
        activeKey = btn.dataset.chartKey;
        draw();
      });
    });

    renderValidationChart(vc);
  };

  draw();
}
