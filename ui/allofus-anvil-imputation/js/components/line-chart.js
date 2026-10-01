/**
 * Line chart (Array / Low-Pass): imputation quality vs. allele frequency, one line per reference panel.
 * Registered as the default validation chart type ("line"). Requires js/components/chart-common.js.
 */
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

VALIDATION_CHART_TYPES.line = { config: lineChartConfig, legendHTML: () => '' };
