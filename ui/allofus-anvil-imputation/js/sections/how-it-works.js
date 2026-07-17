/**
 * Frame 5: how-it-works steps, step connector, and documentation link.
 */
function renderHowItWorksSection(p) {
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
}
