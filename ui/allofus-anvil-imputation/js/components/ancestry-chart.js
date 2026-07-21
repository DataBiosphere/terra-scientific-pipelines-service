/**
 * Renders the ancestry breakdown donut chart as raw SVG paths.
 */
function renderAncestryChart(svgEl, counts, colors) {
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
