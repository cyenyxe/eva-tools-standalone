package uk.ac.ebi.eva.countstats.controller;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.*;
import uk.ac.ebi.eva.countstats.model.CountDto;
import uk.ac.ebi.eva.countstats.service.CountService;

import javax.ws.rs.core.Response;
import java.util.List;

@RestController
@RequestMapping("/v1/countstats")
public class CountController {
    @Autowired
    private CountService countService;

    @GetMapping("status")
    public String welcome() {
        return "Running!!!";
    }

    @PostMapping("count")
    public Response saveCount(@RequestBody CountDto countDto) {
        countService.saveCount(countDto);
        return Response.ok().build();
    }

    @GetMapping("count")
    public List<CountDto> getCount() {
        return countService.getAllCounts();
    }
}
