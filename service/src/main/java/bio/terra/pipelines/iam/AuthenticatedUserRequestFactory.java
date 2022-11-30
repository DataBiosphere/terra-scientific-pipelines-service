package bio.terra.pipelines.iam;

import javax.servlet.http.HttpServletRequest;
import org.springframework.stereotype.Component;

// Making this an interface as I'm not sure what request authentication will look like in mc-terra.
@Component
public interface AuthenticatedUserRequestFactory {

  AuthenticatedUserRequest from(HttpServletRequest servletRequest);
}
