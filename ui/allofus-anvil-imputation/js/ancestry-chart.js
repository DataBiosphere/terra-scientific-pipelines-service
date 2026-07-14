/**
 * Renders the donut chart and tab-specific content from PIPELINES data.
 * Pipeline data lives in js/pipeline-data.js — edit there to change per-pipeline values.
 */

function renderDonutChart(svgEl, counts, colors) {
  while (svgEl.firstChild) svgEl.removeChild(svgEl.firstChild);

  const total = counts.reduce((sum, val) => sum + val, 0);
  const cx = 100, cy = 100, outerR = 70, innerR = 50;
  let angle = 0;

  svgEl.setAttribute('transform', 'rotate(-90)');

  counts.forEach((value, i) => {
    const sweep = (value / total) * 2 * Math.PI;
    const x1 = cx + outerR * Math.cos(angle),    y1 = cy + outerR * Math.sin(angle);
    const x2 = cx + outerR * Math.cos(angle + sweep), y2 = cy + outerR * Math.sin(angle + sweep);
    const x3 = cx + innerR * Math.cos(angle + sweep), y3 = cy + innerR * Math.sin(angle + sweep);
    const x4 = cx + innerR * Math.cos(angle),    y4 = cy + innerR * Math.sin(angle);
    const large = sweep > Math.PI ? 1 : 0;

    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    path.setAttribute('d', [
      `M ${x1} ${y1}`,
      `A ${outerR} ${outerR} 0 ${large} 1 ${x2} ${y2}`,
      `L ${x3} ${y3}`,
      `A ${innerR} ${innerR} 0 ${large} 0 ${x4} ${y4}`,
      `Z`,
    ].join(' '));
    path.setAttribute('fill', colors[i % colors.length]);
    svgEl.appendChild(path);
    angle += sweep;
  });

  const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
  line.setAttribute('x1', '100'); line.setAttribute('y1', '100');
  line.setAttribute('x2', '100'); line.setAttribute('y2', '190');
  line.setAttribute('stroke', '#ADB2BA'); line.setAttribute('stroke-width', '2');
  svgEl.appendChild(line);

  const bulb = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
  bulb.setAttribute('cx', '100'); bulb.setAttribute('cy', '100');
  bulb.setAttribute('r', '3'); bulb.setAttribute('fill', '#ADB2BA');
  svgEl.appendChild(bulb);
}

function renderStepConnector(containerEl, stepCount) {
  const firstY = 110, spacing = 200;
  const lastY = firstY + (stepCount - 1) * spacing;
  const totalH = lastY + 90;

  let circles = '';
  for (let i = 0; i < stepCount; i++) {
    circles += `<circle cx="10" cy="${firstY + i * spacing}" r="5" stroke="#074770" fill="white" stroke-width="2"/>`;
  }

  containerEl.innerHTML = `
    <svg height="${totalH}px" width="20px">
      <line x1="10" x2="10" y1="${firstY}" y2="${lastY}" stroke="#074770" stroke-width="2"/>
      ${circles}
    </svg>`;
}

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
      Scientific Validation
      <div class="validation-subtext">${vc.subtitle}</div>
    </div>
    <div class="validation-chart-wrapper">
      <canvas id="validationChartCanvas"></canvas>
    </div>`;

  if (_validationChart) { _validationChart.destroy(); _validationChart = null; }

  Chart.defaults.font.family = 'Montserrat';

  _validationChart = new Chart(
    document.getElementById('validationChartCanvas').getContext('2d'),
    {
      type: 'line',
      data: {
        labels: vc.labels,
        datasets: vc.datasets.map(ds => ({
          label: ds.label,
          data: ds.data,
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
            title: { display: true, text: vc.xAxisLabel, font: { size: 13, weight: '600' }, color: '#074770', padding: { top: 12 } },
            grid: { color: 'rgba(0,0,0,0.06)' },
            ticks: { font: { size: 13 }, color: '#333F52' },
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

function renderPipeline(pipelineKey) {
  const p = PIPELINES[pipelineKey];
  if (!p) return;

  // Frame 4: reference panel section
  document.getElementById('genome-overview').innerHTML = p.genomeOverviewHTML;
  document.getElementById('total-genomes-count').textContent = p.totalGenomesCount;
  document.getElementById('total-genomes-label').innerHTML = p.totalGenomesLabelHTML;

  document.getElementById('ancestry-list').innerHTML =
    p.ancestryRows.map(r =>
      `<li><span class="ancestry-count">${r.count}</span> ${r.label} <span class="ancestry-percent">${r.percent}</span></li>`
    ).join('') +
    `<li class="ancestry-method">${p.ancestryNoteHTML}</li>`;

  const chartRows = [...p.ancestryRows]
    .filter(r => parseFloat(r.percent) >= 1)
    .sort((a, b) => parseInt(b.count.replace(/,/g, ''), 10) - parseInt(a.count.replace(/,/g, ''), 10));
  const counts = chartRows.map(r => parseInt(r.count.replace(/,/g, ''), 10));
  const colors = chartRows.map(r => r.color);
  renderDonutChart(document.getElementById('donutChart'), counts, colors);

  // Frame 5: how-it-works steps
  const stepsEl = document.getElementById('how-steps');
  stepsEl.innerHTML =
    p.howItWorksSteps.map(s => `
      <div class="how-box">
        <div>
          <h3>${s.title}</h3>
          <div>${s.bodyHTML}</div>
        </div>
        <img src="${s.img}" alt="${s.alt}" />
      </div>`
    ).join('') +
    `<a class="learn-more-btn" href="${p.docsUrl}" target="_blank">
      Documentation
      <svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round" style="flex-shrink:0">
        <path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6"/>
        <polyline points="15 3 21 3 21 9"/>
        <line x1="10" y1="14" x2="21" y2="3"/>
      </svg>
    </a>`;

  renderStepConnector(document.getElementById('step-connector'), p.howItWorksSteps.length);

  // Scientific validation (pipeline-specific, optional)
  renderValidationSection(p);

  // Frame pricing
  document.getElementById('frame-pricing').innerHTML = `
    <div class="pricing-header">
      Pricing
      <div class="pricing-subtext">Simple, transparent per-sample pricing — pay only for what you need.</div>
    </div>
    <div class="pricing-price">
      <span class="pricing-amount">${p.pricePerSample}</span>
      <span class="pricing-unit">per sample</span>
    </div>
    <div class="pricing-disclaimer">
      <strong>Special Discounts Available:</strong> Discounts are available for certain beneficial
      activities. Opportunities to apply for discounts are available within the quota purchasing
      process. For more information, please inquire at
      <a href="mailto:data-science-services-support@broadinstitute.org">data-science-services-support@broadinstitute.org</a>.
    </div>`;
}

function initTabs() {
  const tabsContainer = document.querySelector('.pipeline-tabs');
  const content = document.getElementById('pipeline-content');
  const pipelineKeys = Object.keys(PIPELINES);
  const defaultKey = 'lowpass';
  let currentIndex = pipelineKeys.indexOf(defaultKey);

  // Build tab buttons from pipeline data so descriptions stay in one place
  tabsContainer.innerHTML = pipelineKeys.map((key, i) => {
    const p = PIPELINES[key];
    return `<button class="tab-btn${i === currentIndex ? ' active' : ''}" data-tab="${key}">
      <div class="tab-name">${p.name}</div>
      <div class="tab-desc">${p.tabDescription}</div>
    </button>`;
  }).join('');

  const tabBtns = Array.from(tabsContainer.querySelectorAll('.tab-btn'));

  tabBtns.forEach((btn, newIndex) => {
    btn.addEventListener('click', () => {
      if (newIndex === currentIndex) return;

      const goingRight = newIndex > currentIndex;
      const outClass = goingRight ? 'slide-exit-left' : 'slide-exit-right';
      const inClass  = goingRight ? 'slide-enter-right' : 'slide-enter-left';

      tabBtns.forEach(b => b.classList.remove('active'));
      btn.classList.add('active');

      content.classList.add(outClass);
      content.addEventListener('animationend', () => {
        content.classList.remove(outClass);
        currentIndex = newIndex;
        renderPipeline(btn.dataset.tab);
        content.classList.add(inClass);
        content.addEventListener('animationend', () => {
          content.classList.remove(inClass);
        }, { once: true });
      }, { once: true });
    });
  });

  renderPipeline(defaultKey);

  // Shrink tabs once the intro section scrolls out of view
  const tabsWrapper = document.getElementById('pipeline-tabs-wrapper');
  new IntersectionObserver(
    ([entry]) => tabsWrapper.classList.toggle('is-stuck', !entry.isIntersecting),
    { threshold: 0 }
  ).observe(document.getElementById('product-selection'));
}

initTabs();
