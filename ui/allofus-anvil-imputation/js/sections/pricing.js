/**
 * Pricing: interactive price calculator (academic/nonprofit eligibility + sample count)
 * plus a discount disclaimer.
 */
const IMPUTATION_UI_URL_BASE = 'http://localhost:3000/#pipelines/imputation/run';

function renderPricingSection(p) {
  const container = document.getElementById('frame-pricing');

  container.innerHTML = `
    <div class="pricing-header">
      How much does ${p.name} cost?
    </div>
    <div class="pricing-calculator">
      <div class="pricing-intro">Calculate your price per sample by answering a few questions.</div>
      <div class="pricing-eligibility">
        <div class="pricing-eligibility-label">Do any of the following apply to your work?</div>
        <label class="pricing-checkbox-label">
          <input type="checkbox" class="pricing-eligibility-checkbox">
          I am part of an academic or non-profit organization
        </label>
        <label class="pricing-checkbox-label">
          <input type="checkbox" class="pricing-eligibility-checkbox">
          The work I am doing is for non-profit activities
        </label>
      </div>
      <div class="pricing-sample-count">
        <label for="pricing-sample-count-input">How many samples do you plan to impute?</label>
        <input type="number" id="pricing-sample-count-input" min="1" step="1" placeholder="e.g. 1000">
      </div>
      <button type="button" class="pricing-calculate-btn">See price</button>
      <div class="pricing-result" style="display: none;">
        <span class="pricing-amount"></span>
        <span class="pricing-unit">per sample</span>
      </div>
      <a class="pricing-purchase-btn" href="#" target="_blank" style="display: none;">Purchase</a>
      <div class="pricing-purchase-note" style="display: none;">If you have not yet registered, you will be prompted to create an account. Creating an account takes less than 1 minute.</div>
    </div>
    <div class="pricing-disclaimer">
      <strong>Alternative Pricing Available:</strong> Alternative pricing is available for non-profit activities
      conducted at not-for-profit organizations and for science-at-scale (large-scale purchases). To inquire about
      eligibility for science-at-scale pricing or for other questions, please reach out to
      <a href="mailto:data-science-services-support@broadinstitute.org">data-science-services-support@broadinstitute.org</a>.
    </div>`;

  const [orgCheckbox, activitiesCheckbox] = container.querySelectorAll('.pricing-eligibility-checkbox');
  const sampleCountInput = container.querySelector('#pricing-sample-count-input');
  const calculateBtn = container.querySelector('.pricing-calculate-btn');
  const result = container.querySelector('.pricing-result');
  const resultAmount = result.querySelector('.pricing-amount');
  const purchaseBtn = container.querySelector('.pricing-purchase-btn');
  const purchaseNote = container.querySelector('.pricing-purchase-note');

  calculateBtn.addEventListener('click', () => {
    const sampleCount = parseInt(sampleCountInput.value, 10);
    if (!sampleCount || sampleCount < 1) {
      sampleCountInput.classList.add('input-error');
      result.style.display = 'none';
      purchaseBtn.style.display = 'none';
      purchaseNote.style.display = 'none';
      return;
    }
    sampleCountInput.classList.remove('input-error');

    const nonProfitOrganization = orgCheckbox.checked;
    const nonProfitActivities = activitiesCheckbox.checked;

    // Bulk sample counts don't yet affect the rate, but the count is collected here
    // so pricing tiers can be introduced later without changing this flow.
    const isNonProfit = nonProfitOrganization || nonProfitActivities;
    const price = isNonProfit ? p.priceNonProfit : p.priceForProfit;
    resultAmount.textContent = `$${price.toFixed(2)}`;
    result.style.display = '';

    const purchaseParams = new URLSearchParams({
      nonProfitActivities: String(nonProfitActivities),
      nonProfitOrganization: String(nonProfitOrganization),
      pipeline: p.pipelineKey,
    });
    purchaseBtn.href = `${IMPUTATION_UI_URL_BASE}?${purchaseParams.toString()}`;
    purchaseBtn.style.display = '';
    purchaseNote.style.display = '';

    trackEvent('priceCalculated', {
      pipeline: p.pipelineKey,
      sampleCount: sampleCount,
      isNonprofit: isNonProfit,
      computedPrice: price,
    });
  });

  purchaseBtn.addEventListener('click', () => {
    trackEvent('purchaseClicked');
  });
}
