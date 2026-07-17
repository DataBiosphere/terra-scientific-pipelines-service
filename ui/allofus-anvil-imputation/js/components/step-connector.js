/**
 * Renders the vertical dotted/circle connector line alongside the how-it-works steps.
 */
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
