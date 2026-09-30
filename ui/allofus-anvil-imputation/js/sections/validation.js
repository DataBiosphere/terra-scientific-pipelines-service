/**
 * Scientific validation section (pipeline-specific, optional — hidden when p.validationCharts is absent).
 * Renders one toggle button per entry in p.validationCharts (e.g. SNP / INDEL), plus an optional
 * preprint pill below the chart (p.validationPreprint).
 *
 * Two chart types are supported, selected by vc.chartType:
 *   - (default)   line chart, e.g. imputation quality (r²) vs. allele frequency (Array / Low-Pass)
 *   - "dumbbell"  grouped paired-marker chart, e.g. F1 score of reference-panel vs. imputed calls
 *                 per variant class within each superpopulation, with a shaded per-sample range (SV)
 */
let _validationChart = null;

const VALIDATION_TEXT_COLOR = '#333F52';
const VALIDATION_AXIS_TITLE_COLOR = '#074770';
const VALIDATION_GRID_COLOR = 'rgba(0,0,0,0.06)';

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

function validationAxisTitle(text, padding) {
  return { display: true, text, font: { size: 13, weight: '600' }, color: VALIDATION_AXIS_TITLE_COLOR, padding };
}

function hexToRgba(hex, alpha) {
  const n = parseInt(hex.replace('#', ''), 16);
  return `rgba(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255}, ${alpha})`;
}

/* ---------------------------------------------------------------------------
 * Line chart (Array / Low-Pass): quality vs. allele frequency
 * ------------------------------------------------------------------------- */
function lineChartConfig(vc) {
  // Supports two schemas:
  //   - vc.labels: [...] + plain-number data — one label per point (simple case)
  //   - vc.tickLabels: {position: label} + {x, y} point data — lets datasets carry more
  //     points than there are axis labels (e.g. unlabeled points between labeled ticks,
  //     or extra data points beyond the labeled range)
  const tickLabels = vc.tickLabels || Object.fromEntries(vc.labels.map((lbl, i) => [i + 1, lbl]));
  const tickPositions = Object.keys(tickLabels).map(Number);
  const axisType = vc.xAxisType || 'linear';

  const normalizedDatasets = vc.datasets.map(ds => ({
    label: ds.label,
    data: typeof ds.data[0] === 'object' ? ds.data : ds.data.map((y, i) => ({ x: i + 1, y })),
    borderColor: ds.color,
    backgroundColor: ds.dashed ? 'transparent' : 'rgba(7, 71, 112, 0.06)',
    borderWidth: ds.dashed ? 2 : 3,
    borderDash: ds.dashed ? [6, 4] : [],
    pointRadius: 4,
    pointHoverRadius: 5,
    pointBackgroundColor: 'white',
    pointBorderColor: ds.color,
    pointBorderWidth: 2.5,
    fill: !ds.dashed,
    tension: 0.35,
    clip: false,
  }));

  const allX = normalizedDatasets.flatMap(ds => ds.data.map(pt => pt.x)).concat(tickPositions);
  const dataMin = Math.min(...allX);
  const dataMax = Math.max(...allX);

  return {
    type: 'line',
    data: { datasets: normalizedDatasets },
    options: {
      responsive: true,
      maintainAspectRatio: true,
      events: [],
      plugins: {
        legend: {
          position: 'top',
          labels: { font: { size: 14 }, usePointStyle: true, padding: 24 },
        },
        tooltip: {
          callbacks: {
            label: ctx => ` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)}`,
          },
        },
      },
      scales: {
        x: {
          type: axisType,
          min: dataMin,
          max: dataMax,
          afterBuildTicks: axis => {
            axis.ticks = tickPositions.map(v => ({ value: v }));
          },
          title: validationAxisTitle(vc.xAxisLabel, { top: 12 }),
          grid: { color: VALIDATION_GRID_COLOR },
          ticks: {
            font: { size: 13 },
            color: VALIDATION_TEXT_COLOR,
            callback: value => tickLabels[value] !== undefined ? tickLabels[value] : '',
          },
        },
        y: {
          title: validationAxisTitle(vc.yAxisLabel),
          min: 0, max: 1,
          grid: { color: VALIDATION_GRID_COLOR },
          ticks: { font: { size: 13 }, color: VALIDATION_TEXT_COLOR, stepSize: 0.2 },
        },
      },
    },
  };
}

/* ---------------------------------------------------------------------------
 * Dumbbell chart (SV): paired reference / imputed markers per series within each group
 *
 * Each group (e.g. a superpopulation) sits at integer x = 1..N. Within a group, the series
 * (e.g. SNV / INDEL / SV) are fanned out around the group centre, and each series shows an
 * open marker (reference) and a filled marker (imputed) joined by a stepped connector, with a
 * translucent band behind each marker spanning the per-sample range.
 * ------------------------------------------------------------------------- */
const DUMBBELL_MARKER_COLOR = '#333F52';
const DUMBBELL_MARKER_RADIUS = 5;
const DUMBBELL_SERIES_SPAN = 0.56;     // total x-span the series fan out over within a group
const DUMBBELL_PAIR_HALF_GAP = 0.085;  // half the x-distance between the paired markers
const DUMBBELL_BAND_HALF_WIDTH = 4;    // px — range band is 8px wide
const DUMBBELL_BAND_ALPHA = 0.4;

function dumbbellSeriesOffsets(n) {
  if (n === 1) return [0];
  return Array.from({ length: n }, (_, i) => -DUMBBELL_SERIES_SPAN / 2 + i * DUMBBELL_SERIES_SPAN / (n - 1));
}

function dumbbellXPositions(vc, groupIndex, seriesIndex) {
  const center = groupIndex + 1 + dumbbellSeriesOffsets(vc.series.length)[seriesIndex];
  return { center, reference: center - DUMBBELL_PAIR_HALF_GAP, imputed: center + DUMBBELL_PAIR_HALF_GAP };
}

function drawRoundedRect(ctx, x, y, w, h, r) {
  const radius = Math.min(r, w / 2, h / 2);
  ctx.beginPath();
  ctx.moveTo(x + radius, y);
  ctx.lineTo(x + w - radius, y);
  ctx.arcTo(x + w, y, x + w, y + radius, radius);
  ctx.lineTo(x + w, y + h - radius);
  ctx.arcTo(x + w, y + h, x + w - radius, y + h, radius);
  ctx.lineTo(x + radius, y + h);
  ctx.arcTo(x, y + h, x, y + h - radius, radius);
  ctx.lineTo(x, y + radius);
  ctx.arcTo(x, y, x + radius, y, radius);
  ctx.closePath();
}

function drawDumbbellRangeBand(ctx, px, range, yScale, color) {
  if (!range) return;
  const top = yScale.getPixelForValue(Math.max(range[0], range[1]));
  const bottom = yScale.getPixelForValue(Math.min(range[0], range[1]));
  const height = Math.max(bottom - top, 2 * DUMBBELL_BAND_HALF_WIDTH); // never thinner than it is wide
  const centerY = (top + bottom) / 2;
  ctx.fillStyle = hexToRgba(color, DUMBBELL_BAND_ALPHA);
  drawRoundedRect(ctx, px - DUMBBELL_BAND_HALF_WIDTH, centerY - height / 2, 2 * DUMBBELL_BAND_HALF_WIDTH, height, DUMBBELL_BAND_HALF_WIDTH);
  ctx.fill();
}

// Draws group separators, range bands, and connectors underneath the marker datasets.
function dumbbellDecorationsPlugin(vc) {
  return {
    id: 'dumbbellDecorations',
    beforeDatasetsDraw(chart) {
      const { ctx, chartArea, scales: { x, y } } = chart;
      ctx.save();

      ctx.strokeStyle = VALIDATION_GRID_COLOR;
      ctx.lineWidth = 1;
      for (let i = 1; i < vc.groups.length; i++) {
        const px = x.getPixelForValue(i + 0.5);
        ctx.beginPath(); ctx.moveTo(px, chartArea.top); ctx.lineTo(px, chartArea.bottom); ctx.stroke();
      }

      vc.series.forEach((s, si) => s.data.forEach((pt, gi) => {
        const xs = dumbbellXPositions(vc, gi, si);
        drawDumbbellRangeBand(ctx, x.getPixelForValue(xs.reference), pt.referenceRange, y, s.color);
        drawDumbbellRangeBand(ctx, x.getPixelForValue(xs.imputed), pt.imputedRange, y, s.color);
      }));

      ctx.strokeStyle = DUMBBELL_MARKER_COLOR;
      ctx.lineWidth = 2.5;
      ctx.lineJoin = 'round';
      ctx.lineCap = 'round';
      vc.series.forEach((s, si) => s.data.forEach((pt, gi) => {
        const xs = dumbbellXPositions(vc, gi, si);
        const x0 = x.getPixelForValue(xs.reference);
        const xm = x.getPixelForValue(xs.center);
        const x1 = x.getPixelForValue(xs.imputed);
        const y0 = y.getPixelForValue(pt.reference);
        const y1 = y.getPixelForValue(pt.imputed);
        ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(xm, y0); ctx.lineTo(xm, y1); ctx.lineTo(x1, y1); ctx.stroke();
      }));

      ctx.restore();
    },
  };
}

function dumbbellChartConfig(vc) {
  const groups = vc.groups;
  const markerBase = {
    type: 'scatter',
    pointRadius: DUMBBELL_MARKER_RADIUS,
    pointHoverRadius: DUMBBELL_MARKER_RADIUS,
    pointBorderColor: DUMBBELL_MARKER_COLOR,
    pointBorderWidth: 2,
    showLine: false,
    clip: false,
  };

  const datasets = vc.series.flatMap((s, si) => [
    {
      ...markerBase,
      label: `${s.label} — ${vc.referenceLabel}`,
      pointBackgroundColor: 'white',
      data: s.data.map((pt, gi) => ({ x: dumbbellXPositions(vc, gi, si).reference, y: pt.reference })),
    },
    {
      ...markerBase,
      label: `${s.label} — ${vc.imputedLabel}`,
      pointBackgroundColor: DUMBBELL_MARKER_COLOR,
      data: s.data.map((pt, gi) => ({ x: dumbbellXPositions(vc, gi, si).imputed, y: pt.imputed })),
    },
  ]);

  // Chart.js 4 reads the label colour from each legend item, so every item sets fontColor.
  const legendItemBase = { hidden: false, fontColor: VALIDATION_TEXT_COLOR };
  const legendItems = [
    ...vc.series.map(s => ({
      ...legendItemBase, text: s.label, pointStyle: 'rectRounded',
      fillStyle: hexToRgba(s.color, DUMBBELL_BAND_ALPHA), strokeStyle: hexToRgba(s.color, DUMBBELL_BAND_ALPHA), lineWidth: 0,
    })),
    { ...legendItemBase, text: vc.referenceLabel, pointStyle: 'circle', fillStyle: 'white', strokeStyle: DUMBBELL_MARKER_COLOR, lineWidth: 2 },
    { ...legendItemBase, text: vc.imputedLabel, pointStyle: 'circle', fillStyle: DUMBBELL_MARKER_COLOR, strokeStyle: DUMBBELL_MARKER_COLOR, lineWidth: 2 },
  ];

  const yDecimals = vc.yTickDecimals !== undefined ? vc.yTickDecimals : 2;

  return {
    type: 'scatter',
    data: { datasets },
    plugins: [dumbbellDecorationsPlugin(vc)],
    options: {
      responsive: true,
      maintainAspectRatio: true,
      aspectRatio: 1.75,
      events: [],
      layout: { padding: { top: 4, right: 8 } },
      plugins: {
        legend: {
          position: 'top',
          labels: { font: { size: 14 }, color: VALIDATION_TEXT_COLOR, usePointStyle: true, padding: 20, generateLabels: () => legendItems },
        },
        tooltip: { enabled: false },
      },
      scales: {
        x: {
          type: 'linear',
          min: 0.5,
          max: groups.length + 0.5,
          afterBuildTicks: axis => {
            axis.ticks = groups.map((_, i) => ({ value: i + 1 }));
          },
          title: validationAxisTitle(vc.xAxisLabel, { top: 12 }),
          grid: { display: false },
          border: { color: 'rgba(0,0,0,0.12)' },
          ticks: {
            font: { size: 13, weight: '600' },
            color: VALIDATION_TEXT_COLOR,
            callback: value => groups[Math.round(value) - 1] || '',
          },
        },
        y: {
          min: vc.yMin,
          max: vc.yMax,
          title: validationAxisTitle(vc.yAxisLabel),
          grid: { color: VALIDATION_GRID_COLOR },
          border: { display: false },
          ticks: { font: { size: 13 }, color: VALIDATION_TEXT_COLOR, stepSize: vc.yStepSize, callback: value => Number(value).toFixed(yDecimals) },
        },
      },
    },
  };
}

/* ------------------------------------------------------------------------- */

function renderValidationChart(vc) {
  const canvas = document.getElementById('validationChartCanvas');
  if (_validationChart) { _validationChart.destroy(); _validationChart = null; }

  Chart.defaults.font.family = 'Montserrat';

  const config = vc.chartType === 'dumbbell' ? dumbbellChartConfig(vc) : lineChartConfig(vc);
  _validationChart = new Chart(canvas.getContext('2d'), config);
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
    const wrapperClass = vc.chartType === 'dumbbell' ? 'validation-chart-wrapper validation-chart-wrapper--wide' : 'validation-chart-wrapper';
    container.innerHTML = `
      <div class="validation-header">
        Has the imputation service been scientifically validated?
        <div class="validation-subtext">${vc.subtitle}</div>
      </div>
      ${renderChartToggle()}
      <div class="${wrapperClass}">
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
