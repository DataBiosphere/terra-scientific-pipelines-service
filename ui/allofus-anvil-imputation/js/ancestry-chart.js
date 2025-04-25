
const data = [254416, 101982, 90553, 44627, 13226, 9710];
const colors = ['#2A51B3', '#46A3E9', '#F6BD41', '#5CC88D', '#80C6EC', '#775FE5'];

const total = data.reduce((sum, val) => sum + val, 0);
const centerX = 100;
const centerY = 100;
const outerRadius = 70;
const innerRadius = 50;
let cumulativeAngle = 0;

const svg = document.getElementById('donutChart');
svg.setAttribute('transform', 'rotate(-90)')

data.forEach((value, index) => {
  const sliceAngle = (value / total) * 2 * Math.PI;
  const x1 = centerX + outerRadius * Math.cos(cumulativeAngle);
  const y1 = centerY + outerRadius * Math.sin(cumulativeAngle);
  const x2 = centerX + outerRadius * Math.cos(cumulativeAngle + sliceAngle);
  const y2 = centerY + outerRadius * Math.sin(cumulativeAngle + sliceAngle);

  const x3 = centerX + innerRadius * Math.cos(cumulativeAngle + sliceAngle);
  const y3 = centerY + innerRadius * Math.sin(cumulativeAngle + sliceAngle);
  const x4 = centerX + innerRadius * Math.cos(cumulativeAngle);
  const y4 = centerY + innerRadius * Math.sin(cumulativeAngle);

  const largeArc = sliceAngle > Math.PI ? 1 : 0;

  const pathData = [
    `M ${x1} ${y1}`,
    `A ${outerRadius} ${outerRadius} 0 ${largeArc} 1 ${x2} ${y2}`,
    `L ${x3} ${y3}`,
    `A ${innerRadius} ${innerRadius} 0 ${largeArc} 0 ${x4} ${y4}`,
    `Z`
  ].join(' ');

  const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
  path.setAttribute('d', pathData);
  path.setAttribute('fill', colors[index % colors.length]);
  svg.appendChild(path);

  cumulativeAngle += sliceAngle;
});

// Add line from center to right
const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
line.setAttribute("x1", "100");
line.setAttribute("y1", "100");
line.setAttribute("x2", "100");
line.setAttribute("y2", "190"); // Note SVG is rotated -90 degrees
line.setAttribute("stroke", "#ADB2BA");
line.setAttribute("stroke-width", "2");
svg.appendChild(line);

const bulb = document.createElementNS("http://www.w3.org/2000/svg", "circle");
bulb.setAttribute("cx", "100"); // or "100"
bulb.setAttribute("cy", "100");
bulb.setAttribute("r", "3");
bulb.setAttribute("fill", "#ADB2BA");
svg.appendChild(bulb);
