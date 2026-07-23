/**
 * Minimal Mixpanel event tracking via Bard (https://terra-bard-dev.appspot.com).
 * Visitors here aren't authenticated with Terra, so events are sent unauthenticated
 * with a client-generated, localStorage-persisted anonymous id.
 */
const BARD_ROOT = 'https://terra-bard-dev.appspot.com';
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
