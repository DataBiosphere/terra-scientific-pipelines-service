package bio.terra.pipelines.service;

import org.springframework.stereotype.Service;

@Service
public class EmptyNetworkCallService {

    public void makeANetworkCall() {
        try {
            Thread.sleep(3000);
        } catch (Exception e) {
            // we really, really don't care here
        }
    }
}
