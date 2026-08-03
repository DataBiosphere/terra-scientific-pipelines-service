/**
 * Minimal Mixpanel event tracking via Bard. Points at prod Bard when this page is served
 * from one of its production hostnames, and dev Bard otherwise (dev deployment, local, etc).
 * Visitors here aren't authenticated with Terra, so events are sent unauthenticated
 * with a client-generated, localStorage-persisted anonymous id.
 */
const PROD_HOSTNAMES = [
  'allofus-anvil-imputation.broadinstitute.org',
  'allofus-anvil-imputation.terra.bio',
  'imputation.researchallofus.org',
];
const BARD_ROOT = PROD_HOSTNAMES.includes(window.location.hostname)
  ? 'https://terra-bard-prod.appspot.com'
  : 'https://terra-bard-dev.appspot.com';
const BARD_APP_ID = 'allofus-anvil-imputation-marketing';
const ANON_ID_STORAGE_KEY = 'distinctId';

function getAnonId() {
  let anonId = localStorage.getItem(ANON_ID_STORAGE_KEY);
  if (!anonId) {
    anonId = crypto.randomUUID();
    localStorage.setItem(ANON_ID_STORAGE_KEY, anonId);
  }
  return anonId;
}

function trackEvent(eventName, properties = {}) {
  fetch(`${BARD_ROOT}/api/event`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      event: eventName,
      properties: {
        appId: BARD_APP_ID,
        distinct_id: getAnonId(),
        ...properties
      }
    }),
  }).catch(() => {});
}

trackEvent('pageView', { path: window.location.pathname, referrer: document.referrer });
