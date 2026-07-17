/**
 * Scientific validation section (pipeline-specific, optional — hidden when p.validationChart is absent).
 */
let _validationChart = null;

function renderValidationSection(p) {
  const container = document.getElementById('frame-validation');
  if (!p.validationChart) {
    container.innerHTML = '';
    container.style.display = 'none';
    if (_validationChart) { _validationChart.destroy(); _validationChart = null; }
    return;
  }
  container.style.display = '';

  const vc = p.validationChart;
  container.innerHTML = `
    <div class="validation-header">
      Has the imputation service been scientifically validated?
      <div class="validation-subtext">${vc.subtitle}</div>
    </div>
    <div class="validation-chart-wrapper">
      <canvas id="validationChartCanvas"></canvas>
    </div>`;

  if (_validationChart) { _validationChart.destroy(); _validationChart = null; }

  Chart.defaults.font.family = 'Montserrat';

  // Supports two schemas:
  //   - vc.labels: [...] + plain-number data — one label per point (simple case)
  //   - vc.tickLabels: {position: label} + {x, y} point data — lets datasets carry more
  //     points than there are axis labels (e.g. unlabeled points between labeled ticks)
  const tickLabels = vc.tickLabels || Object.fromEntries(vc.labels.map((lbl, i) => [i + 1, lbl]));
  const tickPositions = Object.keys(tickLabels).map(Number);
  const xMin = Math.min(...tickPositions);
  const xMax = Math.max(...tickPositions);

  _validationChart = new Chart(
    document.getElementById('validationChartCanvas').getContext('2d'),
    {
      type: 'line',
      data: {
        datasets: vc.datasets.map(ds => ({
          label: ds.label,
          data: typeof ds.data[0] === 'object' ? ds.data : ds.data.map((y, i) => ({ x: i + 1, y })),
          borderColor: ds.color,
          backgroundColor: ds.dashed ? 'transparent' : 'rgba(7, 71, 112, 0.06)',
          borderWidth: ds.dashed ? 2 : 3,
          borderDash: ds.dashed ? [6, 4] : [],
          pointRadius: 5,
          pointHoverRadius: 7,
          fill: !ds.dashed,
          tension: 0.35,
        })),
      },
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
            type: 'linear',
            min: xMin,
            max: xMax,
            afterBuildTicks: axis => {
              axis.ticks = tickPositions.map(v => ({ value: v }));
            },
            title: { display: true, text: vc.xAxisLabel, font: { size: 13, weight: '600' }, color: '#074770', padding: { top: 12 } },
            grid: { color: 'rgba(0,0,0,0.06)' },
            ticks: {
              font: { size: 13 },
              color: '#333F52',
              callback: value => tickLabels[value] !== undefined ? tickLabels[value] : '',
            },
          },
          y: {
            title: { display: true, text: vc.yAxisLabel, font: { size: 13, weight: '600' }, color: '#074770' },
            min: 0, max: 1,
            grid: { color: 'rgba(0,0,0,0.06)' },
            ticks: { font: { size: 13 }, color: '#333F52', stepSize: 0.2 },
          },
        },
      },
    }
  );
}
