function renderReferencePanelSection(pipeline) {
  document.getElementById('genome-overview').innerHTML = pipeline.genomeOverviewHTML;
  document.getElementById('total-genomes-count').textContent = pipeline.totalGenomesCount;
  document.getElementById('total-genomes-label').innerHTML = pipeline.totalGenomesLabelHTML;

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
}
