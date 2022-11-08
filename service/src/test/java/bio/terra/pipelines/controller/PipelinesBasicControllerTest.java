package bio.terra.pipelines.controller;

//
// @ContextConfiguration(classes = ExampleController.class)
// @WebMvcTest
// public class ExampleControllerTest {
//  @MockBean ExampleService serviceMock;
//  @MockBean SamUserFactory samUserFactoryMock;
//  @MockBean BearerTokenFactory bearerTokenFactory;
//  @MockBean SamConfiguration samConfiguration;
//  @MockBean SamService samService;
//
//  @Autowired private MockMvc mockMvc;
//
//  private SamUser testUser =
//          new SamUser(
//                  "test@email",
//                  UUID.randomUUID().toString(),
//                  new BearerToken(UUID.randomUUID().toString()));
//
//  @BeforeEach
//  void beforeEach() {
//    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
//  }
//
//  @Test
//  void testGetMessageOk() throws Exception {
//    var example = new Example(testUser.getSubjectId(), "message");
//    when(serviceMock.getExampleForUser(testUser.getSubjectId())).thenReturn(Optional.of(example));
//
//    mockMvc
//            .perform(get("/api/example/v1/message"))
//            .andExpect(status().isOk())
//            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
//            .andExpect(content().string(example.message()));
//  }
//
//  @Test
//  void testGetMessageNotFound() throws Exception {
//    when(serviceMock.getExampleForUser(testUser.getSubjectId())).thenReturn(Optional.empty());
//
//    mockMvc.perform(get("/api/example/v1/message")).andExpect(status().isNotFound());
//  }
// }
