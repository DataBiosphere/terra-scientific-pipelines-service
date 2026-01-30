# Teaspoons email notification templates

This directory contains the email notification templates used by Teaspoons. The templates live in Twilio / SendGrid and 
there is no automated connection between this repository and the Twilio/SendGrid system. Any edits made in that system 
should be copied here. 

To get access to SendGrid message DevOps and ask that they invite `<your-email-username>+dsp-dev-sendgrid@broadinstitute.org` 
and `<your-email-username>+dsp-prod-sendgrid@broadinstitute.org` to Dev and Prod SendGrid respectively.

The HTML files can be previewed by pasting the contents into a preview site such as https://htmledit.squarefree.com.

## Template list
- `job_succeeded.html` - Sent when a job succeeds
- `job_failed.html` - Sent when a job fails
- `user_quota_change_notification.html` - Sent when a user's quota is changed (by calling the Admin endpoint)

## Template Images
Several images and logos are used in the templates. They're included in the `images` directory and referenced in the HTML files. 
To update an image, replace the file in the `images` directory, upload it to the SendGrid media library, and update the image URL in the HTML files to point to the uploaded file.

## Adding/Updating Templates
[This document](https://docs.google.com/document/d/1e2UisAYbW9wXwyI7sWLtYMmX0ILWBP_lvRONBVEbliY/edit?tab=t.0) explains the 
overall process of adding a new notification template and how to test it.

### Updating Existing HTML Templates
For only changing HTML content, you can directly go to the Twilio/SendGrid system and make the changes there.
[This document](https://docs.google.com/document/d/1NSGZoV9Y_wXDfETu4Noqs5PhO7p92XL262ZA3qMStNA/edit?tab=t.0#heading=h.6qwx2xm1aps) 
explains that process and how to test it.

### Notification Parameter Changes
When adding a new notification or changing parameters of an existing notification, changes need to be made to different repos.
More details can be found in [this document](https://docs.google.com/document/d/1e2UisAYbW9wXwyI7sWLtYMmX0ILWBP_lvRONBVEbliY/edit?tab=t.0).

Changes involved:
- terra-helmfile - Add the message ID and its name to Thurloe’s helm configuration in [dev values.yaml](https://github.com/broadinstitute/terra-helmfile/blob/master/charts/thurloe/values.yaml#L33)
  and [prod values.yaml](https://github.com/broadinstitute/terra-helmfile/blob/master/values/app/thurloe/live/prod.yaml)
- workbench Libs - Add the message and its parameters as a new case class in [Notification.scala](https://github.com/broadinstitute/workbench-libs/blob/develop/notifications/src/main/scala/org/broadinstitute/dsde/workbench/model/Notification.scala)
- thurloe - Add the message as a new case in [NotificationMonitor.scala](https://github.com/broadinstitute/thurloe/blob/develop/src/main/scala/thurloe/notification/NotificationMonitor.scala). 
Use the updated workbench-libs version in `build.sbt` to import the message type to Thurloe
- this repo - 
  - Add a corresponding new notification objects at path `service/src/main/java/bio/terra/pipelines/notifications` and use it as needed
  - Copy the HTML content of the template from Twilio/SendGrid and add it to a new file in this repo under `notification-templates/`.
  - Add the images used in the template to the `images` directory if any.
