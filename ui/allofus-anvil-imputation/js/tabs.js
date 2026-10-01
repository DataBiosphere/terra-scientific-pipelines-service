/**
 * Tab switching and per-pipeline render orchestration.
 * Pipeline data lives in js/pipeline-data.js; each section's rendering lives in js/sections/.
 *
 * The selected pipeline is routable via the `?pipeline=` query param, e.g. ?pipeline=sv_imputation.
 * The param accepts a pipeline's `pipelineKey` (matching the Teaspoons UI) or its short PIPELINES key.
 * Adding `&section=pricing` additionally scrolls to that pipeline's pricing calculator on load.
 */
const PIPELINE_QUERY_PARAM = 'pipeline';
const DEFAULT_PIPELINE = 'lowpass';
const SECTION_QUERY_PARAM = 'section';
// Sections that can be linked with `?section=`, mapped to their container element ids.
const ROUTABLE_SECTIONS = { pricing: 'frame-pricing' };

// Resolves the `?pipeline=` value to a PIPELINES key, or null when absent/unrecognised.
function pipelineKeyFromQuery() {
  const value = new URLSearchParams(window.location.search).get(PIPELINE_QUERY_PARAM);
  if (!value) return null;
  return Object.keys(PIPELINES).find(key => key === value || PIPELINES[key].pipelineKey === value) || null;
}

// Resolves the `?section=` value to its container element, or null when absent/unrecognised.
function sectionFromQuery() {
  const value = new URLSearchParams(window.location.search).get(SECTION_QUERY_PARAM);
  return value && ROUTABLE_SECTIONS[value] ? document.getElementById(ROUTABLE_SECTIONS[value]) : null;
}

// Scrolls a section's top to sit just below the sticky tab bar, which is in its compact state once scrolled.
function scrollToSection(section) {
  if (!section || section.offsetParent === null) return; // hidden, e.g. for a Coming Soon pipeline
  const tabsWrapper = document.getElementById('pipeline-tabs-wrapper');
  // Measure the compact bar height with its padding transition suppressed, otherwise the
  // mid-transition height is read and the section lands too low.
  tabsWrapper.classList.add('is-stuck', 'no-transition');
  const top = section.getBoundingClientRect().top + window.scrollY - tabsWrapper.offsetHeight;
  tabsWrapper.classList.remove('no-transition');
  window.scrollTo({ top, behavior: 'auto' });
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

  // Link straight to a section of the selected pipeline, e.g. ?pipeline=sv_imputation&section=pricing
  const linkedSection = sectionFromQuery();
  if (linkedSection) {
    scrollToSection(linkedSection);
    // Re-align once images and fonts have loaded, since they can shift the layout above the section
    if (document.readyState !== 'complete') {
      window.addEventListener('load', () => scrollToSection(linkedSection), { once: true });
    }
  }
}

initTabs();
