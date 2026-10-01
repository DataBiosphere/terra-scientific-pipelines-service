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
 * open marker (reference) and a filled marker (imputed) joined by a stepped connector. When
 * per-sample values are given, each side gets a half-violin (KDE) on the pair's centre line plus
 * small dots for the individual samples; otherwise a translucent band spans the given range.
 * ------------------------------------------------------------------------- */
const DUMBBELL_MARKER_COLOR = '#333F52';
const DUMBBELL_MARKER_RADIUS = 5;
const DUMBBELL_SERIES_SPAN = 0.56;     // total x-span the series fan out over within a group
const DUMBBELL_PAIR_HALF_GAP = 0.05;   // half the x-distance between the paired markers
const DUMBBELL_BAND_HALF_WIDTH = 4;    // px — range band is 8px wide
const DUMBBELL_BAND_ALPHA = 0.4;
const DUMBBELL_GROUP_SHADE = 'rgba(7, 71, 112, 0.045)'; // alternating group background
const DUMBBELL_VIOLIN_MAX_HALF_WIDTH = 16; // px — every violin's widest point reaches this far from the centre line
const DUMBBELL_VIOLIN_OUTLINE_ALPHA = 0.75;
const DUMBBELL_VIOLIN_GRID_POINTS = 40;
const DUMBBELL_SAMPLE_DOT_COLOR = 'rgba(51, 63, 82, 0.45)'; // individual per-sample dots
const DUMBBELL_SAMPLE_DOT_RADIUS = 2.25;
const DUMBBELL_SAMPLE_DOT_JITTER = 2.5;   // px — alternate dots left/centre/right so stacked values stay visible

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

function hasSamples(samples) {
  return Array.isArray(samples) && samples.length >= 2;
}

function sampleRange(samples) {
  return [Math.min(...samples), Math.max(...samples)];
}

// Gaussian KDE with Scott's bandwidth (what seaborn/scipy use by default), evaluated on an even grid
// between the sample min and max (i.e. cut = 0). Returns {values, densities} with densities scaled to max 1.
function kernelDensity(samples) {
  const n = samples.length;
  const mean = samples.reduce((a, b) => a + b, 0) / n;
  const variance = samples.reduce((a, b) => a + (b - mean) ** 2, 0) / (n - 1);
  const bandwidth = Math.sqrt(variance) * Math.pow(n, -0.2);
  const [lo, hi] = sampleRange(samples);
  if (!(bandwidth > 0) || hi - lo === 0) return null; // all samples identical: no spread to draw
  const values = Array.from({ length: DUMBBELL_VIOLIN_GRID_POINTS }, (_, i) => lo + (hi - lo) * i / (DUMBBELL_VIOLIN_GRID_POINTS - 1));
  const densities = values.map(v => samples.reduce((a, s) => a + Math.exp(-0.5 * ((v - s) / bandwidth) ** 2), 0));
  const peak = Math.max(...densities);
  return { values, densities: densities.map(dn => dn / peak) };
}

// Draws one half of a split violin: a flat edge on the pair's centre line, bulging towards `direction`
// (-1 left for the reference side, +1 right for the imputed side).
function drawDumbbellHalfViolin(ctx, centerPx, samples, direction, yScale, color) {
  const kde = kernelDensity(samples);
  ctx.fillStyle = hexToRgba(color, DUMBBELL_BAND_ALPHA);
  ctx.strokeStyle = hexToRgba(color, DUMBBELL_VIOLIN_OUTLINE_ALPHA);
  ctx.lineWidth = 1;
  ctx.lineJoin = 'round';
  if (!kde) {
    // Degenerate distribution: a short flat sliver at the shared value so the side still reads as present
    const yPx = yScale.getPixelForValue(samples[0]);
    ctx.fillRect(Math.min(centerPx, centerPx + direction * DUMBBELL_VIOLIN_MAX_HALF_WIDTH * 0.5), yPx - 1.5, DUMBBELL_VIOLIN_MAX_HALF_WIDTH * 0.5, 3);
    return;
  }
  ctx.beginPath();
  ctx.moveTo(centerPx, yScale.getPixelForValue(kde.values[kde.values.length - 1]));
  for (let i = kde.values.length - 1; i >= 0; i--) {
    ctx.lineTo(centerPx + direction * kde.densities[i] * DUMBBELL_VIOLIN_MAX_HALF_WIDTH, yScale.getPixelForValue(kde.values[i]));
  }
  ctx.lineTo(centerPx, yScale.getPixelForValue(kde.values[0]));
  ctx.closePath();
  ctx.fill();
  ctx.stroke();
}

function drawDumbbellSampleDots(ctx, px, samples, yScale) {
  if (!hasSamples(samples)) return;
  ctx.fillStyle = DUMBBELL_SAMPLE_DOT_COLOR;
  samples.forEach((value, i) => {
    const dx = ((i % 3) - 1) * DUMBBELL_SAMPLE_DOT_JITTER;
    ctx.beginPath();
    ctx.arc(px + dx, yScale.getPixelForValue(value), DUMBBELL_SAMPLE_DOT_RADIUS, 0, Math.PI * 2);
    ctx.fill();
  });
}

// Draws group separators, violins or range bands, sample dots, and connectors underneath the marker datasets.
function dumbbellDecorationsPlugin(vc) {
  return {
    id: 'dumbbellDecorations',
    beforeDatasetsDraw(chart) {
      const { ctx, chartArea, scales: { x, y } } = chart;
      ctx.save();

      // Alternate a faint shade behind every other group so the groups read as separate columns
      ctx.fillStyle = DUMBBELL_GROUP_SHADE;
      for (let i = 1; i < vc.groups.length; i += 2) {
        const left = x.getPixelForValue(i + 0.5);
        const right = x.getPixelForValue(i + 1.5);
        ctx.fillRect(left, chartArea.top, right - left, chartArea.bottom - chartArea.top);
      }

      ctx.strokeStyle = VALIDATION_GRID_COLOR;
      ctx.lineWidth = 1;
      for (let i = 1; i < vc.groups.length; i++) {
        const px = x.getPixelForValue(i + 0.5);
        ctx.beginPath(); ctx.moveTo(px, chartArea.top); ctx.lineTo(px, chartArea.bottom); ctx.stroke();
      }

      vc.series.forEach((s, si) => s.data.forEach((pt, gi) => {
        const xs = dumbbellXPositions(vc, gi, si);
        const centerPx = x.getPixelForValue(xs.center);
        if (hasSamples(pt.referenceSamples)) drawDumbbellHalfViolin(ctx, centerPx, pt.referenceSamples, -1, y, s.color);
        else drawDumbbellRangeBand(ctx, x.getPixelForValue(xs.reference), pt.referenceRange, y, s.color);
        if (hasSamples(pt.imputedSamples)) drawDumbbellHalfViolin(ctx, centerPx, pt.imputedSamples, 1, y, s.color);
        else drawDumbbellRangeBand(ctx, x.getPixelForValue(xs.imputed), pt.imputedRange, y, s.color);
      }));

      vc.series.forEach((s, si) => s.data.forEach((pt, gi) => {
        const xs = dumbbellXPositions(vc, gi, si);
        drawDumbbellSampleDots(ctx, x.getPixelForValue(xs.reference), pt.referenceSamples, y);
        drawDumbbellSampleDots(ctx, x.getPixelForValue(xs.imputed), pt.imputedSamples, y);
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

// Two-row legend rendered above the canvas: marker meaning first, then the series colours.
function dumbbellLegendHTML(vc) {
  const item = (swatch, text) => `<span class="validation-legend-item">${swatch}<span>${text}</span></span>`;
  const markerRow = [
    item('<span class="validation-legend-marker" aria-hidden="true"></span>', vc.referenceLabel),
    item('<span class="validation-legend-marker validation-legend-marker--filled" aria-hidden="true"></span>', vc.imputedLabel),
  ];
  const anySamples = vc.series.some(s => s.data.some(pt => hasSamples(pt.referenceSamples) || hasSamples(pt.imputedSamples)));
  if (anySamples) markerRow.push(item('<span class="validation-legend-dot" aria-hidden="true"></span>', vc.samplesLabel || 'Individual sample'));
  const seriesRow = vc.series.map(s =>
    item(`<span class="validation-legend-swatch" style="background:${hexToRgba(s.color, DUMBBELL_BAND_ALPHA)}" aria-hidden="true"></span>`, s.label)
  );
  return `<div class="validation-legend">
      <div class="validation-legend-row">${markerRow.join('')}</div>
      <div class="validation-legend-row">${seriesRow.join('')}</div>
    </div>`;
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
        legend: { display: false }, // rendered as HTML by dumbbellLegendHTML() so it can span two rows
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
    container.innerHTML = `
      <div class="validation-header">
        Has the imputation service been scientifically validated?
        <div class="validation-subtext">${vc.subtitle}</div>
      </div>
      ${renderChartToggle()}
      <div class="validation-chart-wrapper">
        ${vc.chartType === 'dumbbell' ? dumbbellLegendHTML(vc) : ''}
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
