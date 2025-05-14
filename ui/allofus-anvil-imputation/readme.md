# All of US + AnVIL Imputation Service marketing UI
This code is for user education and landing zones in the All of US + AnVIL Imputation Service.

The UI is plain HTML, CSS, images, and JavaScript that is simple and easy to update and deploy.

The main page is index.html, which maps to https://allofus-anvil-imputation.terra.bio.  A secondary page is acknowledgments.html, which maps to https://terra.bio/allofus-anvil-imputation/acknowledgments/.

Files are hosted in Google Cloud Storage buckets, and configured for HTTPS with a Google-managed SSL certificate via a load balancer.  These resources are configured through Terraform; see [terraform-ap-deployments/teaspoons-imputation-marketing](https://github.com/broadinstitute/terraform-ap-deployments/tree/master/teaspoons-imputation-marketing).

### Development
To develop and test, you can use e.g. [Simple Web Server](https://simplewebserver.org/) to see your local web code rendered as web pages.

### Deployment
To make this UI available at a public URL, use the deployment script like so:

```
cd ui/allofus-anvil-imputation
gcloud auth login # Sign in with your firecloud.org account if deploying to production
python3 deploy.py --environment dev # use "--environment prod"  if deploying to production
```

That will copy files into the appropriate bucket, while ignoring extraneous files, refining cache settings, and adding error handling conveniences.  You can see your deployment at:

* Development: https://allofus-anvil-imputation.dsde-dev.broadinstitute.org
* Production: https://allofus-anvil-imputation.terra.bio

Updates may take up to 10 minutes to appear at the default URL, due to caching.  To see the uncached page, append a random URL parameter to the URL, e.g. https://allofus-anvil-imputation.terra.bio?foo=bar.  To update cache TTL for subsequent deployments, change `max-age=600` in deploy.py.

Please send any questions to scientific-services-support@broadinstitute.org or #imputation-service in Broad Slack.
