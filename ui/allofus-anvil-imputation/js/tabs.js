/**
 * Tab switching and per-pipeline render orchestration.
 * Pipeline data lives in js/pipeline-data.js; each section's rendering lives in js/sections/.
 *
 * The selected pipeline is routable via the `?pipeline=` query param, e.g. ?pipeline=sv_imputation.
 * The param accepts a pipeline's `pipelineKey` (matching the Teaspoons UI) or its short PIPELINES key.
 */
const PIPELINE_QUERY_PARAM = 'pipeline';
const DEFAULT_PIPELINE = 'lowpass';

// Resolves the `?pipeline=` value to a PIPELINES key, or null when absent/unrecognised.
function pipelineKeyFromQuery() {
  const value = new URLSearchParams(window.location.search).get(PIPELINE_QUERY_PARAM);
  if (!value) return null;
  return Object.keys(PIPELINES).find(key => key === value || PIPELINES[key].pipelineKey === value) || null;
}

// Reflects the selected pipeline in the URL without adding a history entry.
function syncPipelineQueryParam(key) {
  const url = new URL(window.location.href);
  url.searchParams.set(PIPELINE_QUERY_PARAM, PIPELINES[key].pipelineKey);
  history.replaceState(null, '', url);
}

function renderPipeline(pipelineKey) {
  const p = PIPELINES[pipelineKey];
  if (!p) return;

  const normalSections = [
    document.getElementById('frame-4'),
    document.getElementById('frame-5'),
    document.getElementById('frame-pricing'),
  ];

  // If a pipeline is Coming Soon, we'll only render the Coming Soon details
  if (p.comingSoon) {
    normalSections.forEach(el => { el.style.display = 'none'; });
    renderValidationSection(p); // not included in Coming Soon pipelines but needed to cleanup between tab switches
    renderComingSoonSection(p);
    return;
  }

  normalSections.forEach(el => { el.style.display = ''; });
  renderComingSoonSection(p);
  renderReferencePanelSection(p);
  renderHowItWorksSection(p);
  renderValidationSection(p);
  renderPricingSection(p);
}

function initTabs() {
  const tabsContainer = document.querySelector('.pipeline-tabs');
  const content = document.getElementById('pipeline-content');
  const pipelineKeys = Object.keys(PIPELINES);
  const defaultKey = pipelineKeyFromQuery() || DEFAULT_PIPELINE;
  let currentIndex = pipelineKeys.indexOf(defaultKey);

  // Build tab buttons from pipeline data so descriptions stay in one place
  tabsContainer.innerHTML = pipelineKeys.map((key, i) => {
    const p = PIPELINES[key];
    return `<button class="tab-btn${i === currentIndex ? ' active' : ''}" data-tab="${key}">
      <div class="tab-name">${p.name}${p.comingSoon ? ' <span class="tab-coming-soon-badge">Coming Soon</span>' : ''}</div>
      <div class="tab-desc">${p.tabDescription}</div>
    </button>`;
  }).join('');

  const tabBtns = Array.from(tabsContainer.querySelectorAll('.tab-btn'));

  tabBtns.forEach((btn, newIndex) => {
    btn.addEventListener('click', () => {
      if (newIndex === currentIndex) return;

      trackEvent('tabSelected', { pipeline: PIPELINES[btn.dataset.tab].pipelineKey });
      syncPipelineQueryParam(btn.dataset.tab);

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
