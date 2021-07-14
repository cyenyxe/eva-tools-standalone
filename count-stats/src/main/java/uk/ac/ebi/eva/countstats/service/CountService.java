package uk.ac.ebi.eva.countstats.service;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import uk.ac.ebi.eva.countstats.model.CountDto;
import uk.ac.ebi.eva.countstats.repository.CountRepository;

import java.util.List;

@Service
public class CountService {
    @Autowired
    private CountRepository countRepository;

    public int saveCount(CountDto countDto) {
        return countRepository.saveCount(countDto);
    }

    public List<CountDto> getAllCounts() {
        return countRepository.getAllCounts();
    }

}
