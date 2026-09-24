const cleanupPlugin = function(system) {
  return {
    components: {
      // we don't need to load different specs here...
      Topbar: () => null,
      // since everything's application/json, the content-type dropdown just takes up space
      contentType: ({ value, contentTypes }) => {
        return system.React.createElement('div', null, value || contentTypes.first())
      }
    }
  }
}

window.onload = function() {
  const swaggerUiEl = document.getElementById('swagger-ui')
  const clientId = swaggerUiEl.dataset.clientId

  // Begin Swagger UI call region
  const ui = SwaggerUIBundle({
    url: 'openapi.yml',
    dom_id: '#swagger-ui',
    deepLinking: true,
    presets: [
      SwaggerUIBundle.presets.apis,
      SwaggerUIStandalonePreset
    ],
    plugins: [
      SwaggerUIBundle.plugins.DownloadUrl,
      cleanupPlugin
    ],
    layout: 'StandaloneLayout',
    defaultModelsExpandDepth: -1, // hide the huge list of schemas under the routes
    defaultModelRendering: 'model', // schema has the descriptions, unlike the example value
    defaultModelExpandDepth: 2, // affects the schema shown for a request or response
    oauth2RedirectUrl: `${window.location.protocol}//${window.location.host}/webjars/swagger-ui-dist/oauth2-redirect.html`
  })
  // End Swagger UI call region

  ui.initOAuth({
    clientId: clientId,
    scopes: "openid",
    additionalQueryStringParams: {prompt: "login"},
    usePkceWithAuthorizationCodeGrant: true
  })

  window.ui = ui
}
