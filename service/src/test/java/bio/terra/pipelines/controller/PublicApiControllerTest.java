package bio.terra.pipelines.controller;

import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.jsonPath;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.pipelines.app.controller.PublicApiController;
import bio.terra.pipelines.generated.model.ApiVersionProperties;
import bio.terra.pipelines.service.StatusService;
import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit.jupiter.SpringExtension;
import org.springframework.test.web.servlet.MockMvc;

@ExtendWith(SpringExtension.class)
@ContextConfiguration(classes = PublicApiController.class)
@WebMvcTest
class PublicApiControllerTest extends BaseTest {

  @Autowired private MockMvc mockMvc;

  @MockBean private ApiVersionProperties versionProperties;
  @MockBean private StatusService statusService;

  @Test
  void testStatus() throws Exception {
    when(statusService.getCurrentStatus()).thenReturn(true);
    this.mockMvc.perform(get("/status")).andExpect(status().isOk());
  }

  @Test
  void testStatusNotOk() throws Exception {
    when(statusService.getCurrentStatus()).thenReturn(false);
    this.mockMvc.perform(get("/status")).andExpect(status().isServiceUnavailable());
  }

  @Test
  void testVersion() throws Exception {
    String gitTag = "0.1.0";
    String gitHash = "abc1234";
    String github = "https://github.com/DataBiosphere/terra-policy-service/tree/0.0.0";
    String build = "0.1.0";

    when(versionProperties.getGitTag()).thenReturn(gitTag);
    when(versionProperties.getGitHash()).thenReturn(gitHash);
    when(versionProperties.getGithub()).thenReturn(github);
    when(versionProperties.getBuild()).thenReturn(build);

    this.mockMvc
        .perform(get("/version"))
        .andExpect(status().isOk())
        .andExpect(jsonPath("$.gitTag").value(gitTag))
        .andExpect(jsonPath("$.gitHash").value(gitHash))
        .andExpect(jsonPath("$.github").value(github))
        .andExpect(jsonPath("$.build").value(build));
  }
}
