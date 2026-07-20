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
    <div class="pricing-disclaimer">
      <strong>Alternative Pricing Available:</strong> Alternative pricing is available for non-profit activities 
      conducted at not-for-profit organizations and for science-at-scale (large-scale purchases). An opportunity 
      to indicate eligibility for alternative pricing for non-profit activities is available within the quota purchasing 
      process. To inquire about eligibility for science-at-scale pricing or for other questions, please reach out to
      <a href="mailto:data-science-services-support@broadinstitute.org">data-science-services-support@broadinstitute.org</a>.
    </div>`;
}
