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

// crypto.randomUUID is only available in secure contexts; fall back to a getRandomValues-based v4 UUID.
function randomUUID() {
  if (typeof crypto !== 'undefined' && typeof crypto.randomUUID === 'function') return crypto.randomUUID();
  const bytes = new Uint8Array(16);
  crypto.getRandomValues(bytes);
  bytes[6] = (bytes[6] & 0x0f) | 0x40;
  bytes[8] = (bytes[8] & 0x3f) | 0x80;
  const hex = Array.from(bytes, b => b.toString(16).padStart(2, '0')).join('');
  return `${hex.slice(0, 8)}-${hex.slice(8, 12)}-${hex.slice(12, 16)}-${hex.slice(16, 20)}-${hex.slice(20)}`;
}

// localStorage can throw (e.g. private browsing, blocked storage); fall back to a per-page-load id.
let _sessionAnonId = null;
function getAnonId() {
  try {
    let anonId = localStorage.getItem(ANON_ID_STORAGE_KEY);
    if (!anonId) {
      anonId = randomUUID();
      localStorage.setItem(ANON_ID_STORAGE_KEY, anonId);
    }
    return anonId;
  } catch (e) {
    if (!_sessionAnonId) _sessionAnonId = randomUUID();
    return _sessionAnonId;
  }
}

// Analytics must never break the page: any failure here is swallowed.
function trackEvent(eventName, properties = {}) {
  try {
    fetch(`${BARD_ROOT}/api/event`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        event: `imputationMarketing:${eventName}`,
        properties: {
          appId: BARD_APP_ID,
          distinct_id: getAnonId(),
          ...properties
        }
      }),
    }).catch(() => {});
  } catch (e) {
    // ignore
  }
}

trackEvent('pageView', { path: window.location.pathname, referrer: document.referrer });
