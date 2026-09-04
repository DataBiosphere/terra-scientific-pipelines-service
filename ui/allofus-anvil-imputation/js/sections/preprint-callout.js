const PREPRINT = {
  url: 'https://www.medrxiv.org/content/10.64898/2026.08.25.26361247v1',
  badge: 'New preprint',
  shortTitle: 'Our 515,579-genome reference panel improves rare-variant imputation',
  linkText: 'Read it on medRxiv',
};

const ARROW_SVG = `<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round" style="flex-shrink:0" aria-hidden="true">
  <line x1="5" y1="12" x2="19" y2="12"/><polyline points="12 5 19 12 12 19"/></svg>`;

function renderPreprintCallout() {
  const slot = document.getElementById('preprint-slot-pill');
  if (!slot) return;

  slot.innerHTML = `<a class="preprint-callout" href="${PREPRINT.url}" target="_blank" rel="noopener">
    <span class="preprint-badge">${PREPRINT.badge}</span>
    <span class="preprint-text">${PREPRINT.shortTitle}</span>
    <span class="preprint-link">${PREPRINT.linkText} ${ARROW_SVG}</span>
  </a>`;
}

renderPreprintCallout();
