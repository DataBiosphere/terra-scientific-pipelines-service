package bio.terra.pipelines;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.web.servlet.config.annotation.CorsRegistry;
import org.springframework.web.servlet.config.annotation.WebMvcConfigurer;

@Configuration
public class WebConfig {

  @Bean
  public WebMvcConfigurer corsConfigurer() {
    return new WebMvcConfigurer() {
      @Override
      public void addCorsMappings(CorsRegistry registry) {
        registry.addMapping("/docs").allowedOrigins("*").allowedMethods("GET", "OPTIONS");
        registry.addMapping("/status").allowedOrigins("*").allowedMethods("GET", "OPTIONS");
        registry.addMapping("/version").allowedOrigins("*").allowedMethods("GET", "OPTIONS");
      }
    };
  }
}
