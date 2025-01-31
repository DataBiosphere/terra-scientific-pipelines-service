# Teaspoons email notification templates

This directory contains the email notification templates used by Teaspoons. The templates live in Twilio / SendGrid and there is no automated connection between this repository and the Twilio/SendGrid system. Any edits made in that system should be copied here.

The HTML files can be previewed by pasting the contents into a preview site such as https://htmledit.squarefree.com.

## Template list
- `job_succeeded.html` - Sent when a job succeeds
- `job_failed.html` - Sent when a job fails


## Updating Templates
Doc explaining the process for changes and how to test - https://docs.google.com/document/d/1e2UisAYbW9wXwyI7sWLtYMmX0ILWBP_lvRONBVEbliY/edit?tab=t.0

### HTML Changes
For these you only have to go to the Twilio/SendGrid system and make the changes there.

### Notification Parameter Changes
Changes need to be made to different repos (in this order):
- Workbench Libs - update https://github.com/broadinstitute/workbench-libs/blob/develop/notifications/src/main/scala/org/broadinstitute/dsde/workbench/model/Notification.scala#L209
- Thurloe - update to version of Workbench Libs in build.sbt generated from above changes and update https://github.com/broadinstitute/thurloe/blob/develop/src/main/scala/thurloe/notification/NotificationMonitor.scala#L454
- This repo - update the notification objects in this repo.
