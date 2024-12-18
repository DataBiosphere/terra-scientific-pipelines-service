package bio.terra.pipelines.notifications;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.springframework.beans.factory.annotation.Autowired;

class PubsubServiceTest extends BaseEmbeddedDbTest {
  @Autowired PubsubService pubsubService;
}
