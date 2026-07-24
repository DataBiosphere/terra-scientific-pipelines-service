function renderComingSoonSection(pipeline) {
  const container = document.getElementById('frame-coming-soon');
  const notice = pipeline.comingSoon;
  if (!notice) {
    container.innerHTML = '';
    container.style.display = 'none';
    return;
  }
  container.style.display = '';

  container.innerHTML = `
    <img class="coming-soon-img" src="img/dna-orbit.png" alt="" />
    <div class="coming-soon-message">${notice.message}</div>
    <a class="learn-more-btn" href="${notice.signupUrl}" target="_blank">${notice.signupLabel}</a>`;
}
