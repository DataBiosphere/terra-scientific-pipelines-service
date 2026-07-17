/**
 * Frame pricing: per-sample price, quota note, and discount disclaimer.
 */
function renderPricingSection(p) {
  document.getElementById('frame-pricing').innerHTML = `
    <div class="pricing-header">
      How much does ${p.name} cost?
    </div>
    <div class="pricing-price">
      <span class="pricing-amount">${p.pricePerSample}</span>
      <span class="pricing-unit">per sample</span>
    </div>
    <div class="pricing-quota-note">
      Quota expires one year after issuance. For multi-year studies planning to impute in multiple years,
      we recommend buying what you need today — you can always easily purchase more via credit card.
    </div>
    <div class="pricing-disclaimer">
      <strong>Special Discounts Available:</strong> Discounts are available for certain beneficial
      activities. An opportunity to indicate eligibility for discounts is available within the quota purchasing
      process. For more information, please inquire at
      <a href="mailto:data-science-services-support@broadinstitute.org">data-science-services-support@broadinstitute.org</a>.
    </div>`;
}
