/**
 * Handle checking the status of a job and
 * updating the progress bar on the waiting page.
 */
async function handleJobStatus(jobStatusURL, sessionId) {
  let checks = 0;
  let jobStatus = "PENDING";
  let redirectUrl = "";

  await new Promise((resolve) => setTimeout(resolve, 1000));

  // Pending loop
  while (jobStatus == "PENDING") {
    let response = await fetch(jobStatusURL);
    var data = await response.json();
    jobStatus = data.status;
    if (jobStatus == "PENDING") {
      // wait 10 seconds
      await new Promise((resolve) => setTimeout(resolve, 10000));
    }
  }

  // Running loop
  let last_stage_index = 0;
  while (jobStatus == "RUNNING") {
    checks += 1;
    let response = await fetch(jobStatusURL);
    var data = await response.json();
    jobStatus = data.status;
    redirectUrl = data.redirect_url ? data.redirect_url : "";
    if (jobStatus == "RUNNING") {
      let stage_index = data.stage_index + 1;
      let stage_count = data.stage_count;
      let gap = (1 / stage_count) * 100;
      if (stage_index !== last_stage_index) {
        checks = 0;
        last_stage_index = stage_index;
      }
      let percent = ((stage_index-0.5) / stage_count) * 100 + (1 - (1/(1+(checks*0.25)))) * gap;
      document.getElementById("progress-bar").style.width = percent + "%";
      document.getElementById("progress-bar").ariaValueNow = percent;

      document.getElementById("job-status-text").innerHTML = `<i>Your job is currently running. Stage ${stage_index} of ${stage_count}.</i>`;
      // wait 10 seconds
      await new Promise((resolve) => setTimeout(resolve, 10000));
    }
  };
  // Remove spinner
  document.getElementById("loading-spinner").setAttribute("style", "display:none !important");

  if (jobStatus == "FAILURE") {
    document.getElementById("progress-bar").style.width = "100%";
    document.getElementById("progress-bar").ariaValueNow = 100;
    document.getElementById("progress-bar").classList.add("bg-danger");
    document.getElementById("error-section").style.display = "block";
    document.getElementById("job-status-text").innerHTML = `<i>Your submission has failed. Please see the error message below.</i>`;
    handleError(data, sessionId);
  } else if (jobStatus == "SUCCESS") {
    document.getElementById("progress-bar").style.width = "100%";
    document.getElementById("progress-bar").ariaValueNow = 100;
    document.getElementById("progress-bar").classList.add("bg-success");
    document.getElementById("job-status-text").innerHTML = `<i>Your results are ready!</i>`;

    document.getElementById("success-section").style.display = "block";
    document.getElementById("success-button").href = redirectUrl;
  }
}

function handleError(data, sessionId) {
  document.getElementById("error-title").innerHTML = data.error_title;
  document.getElementById("error-message").innerHTML = `<code>${data.error_message}</code>`;
  
  if (data.status_code >= 400 && data.status_code < 500) {
    document.getElementById("error-subtitle").innerHTML = "This error is due to an issue with your form input or your uploaded files. Please see the error message below for more details.";
  } else if (data.status_code >= 500) {
    document.getElementById("error-subtitle").innerHTML = "This error is due to an issue with the server. Please see the error message below for more details, and report this to our system administrator using the link below.";
    document.getElementById("error-contact-section").style.display = "block";
    document.getElementById("error-contact-link").href = buildMailToLink(data, sessionId);
  }

  if (data.payload) {
    document.getElementById("error-payload").innerHTML = `<pre><code>${JSON.stringify(data.payload, null, 2)}</code></pre>`;
  }
}

function buildMailToLink(data, sessionId) {
  const subject = `LocusFocus Error Report [${sessionId}]`;
  const body = `Hello,
    I received a server error from LocusFocus from submitting a job. Please see the details below.

    Details:
    - Session ID: ${sessionId}
    - Error Title: ${data.error_title}
    - Error Message: ${data.error_message}
    - Payload: ${data.payload}

    Thanks!
  `;
  const link = `mailto:mackenzie.frew@sickkids.ca?subject=${encodeURIComponent(subject)}&body=${encodeURIComponent(body)}`;
  return link;
}